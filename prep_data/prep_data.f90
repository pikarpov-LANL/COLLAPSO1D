program read 
    !*********************************************************              
    !                                                        *              
    !  This program reads data from models constructed by    *              
    !  Stan Woosley of Supernova precursors.                 *              
    !  It then sets up a file that can be read by our        *              
    !  1d supernova model.                                   *              
    !                                                        *              
    !*********************************************************              
    !                                                                             
          implicit double precision (a-h, o-z) 
    !                                                                       
          parameter(utemp=1e9) 
          parameter(udens=2e6) 
          parameter(uvel=1e8) 
          parameter(udist=1e9) 
          parameter(utime=1e1) 
          parameter(uergg=1e16) 
          parameter(pi43=3.14159*4.0/3.0) 
    !                                                                       
          parameter(idim=4000) 
    !                                                                       
          common /celle/ x(0:idim),v(0:idim) 
          common /cellc/ u(idim),rho(idim),ye(idim),q(idim),dq(idim) 
          common /numb/ ncell 
          common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
         &               unue(idim),unueb(idim),unux(idim)                  
          common /eosq / pr(idim), vsound(idim), u2(idim), vsmax 
          common /carac/ deltam(idim), abar(idim) 
          common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
          common /freez/ ufreez(idim) 
          common /tempe/ temp(idim) 
          common /cent/ dj(idim) 
          common /abun/ ycc 
    
          double precision yccin(idim,20),rnorm 
          real ycc(idim,20) 
    
          double precision rhocgs, tkelv, yej, abarj, ucgs, pcgs, xpj, xnj 
          double precision f3, enclmass(0:idim) 
          double precision mcut 
          double precision, allocatable :: vel(:),rad(:),dens(:),t9(:),     &
         &             yel(:),ab(:),omega(:),press(:)                       
          double precision maxrad, deltam_growth, enclmass_conv_cutoff 
          integer max,j,i,izone,ieos
          integer conv_grid,conv_grid_goal,grid_goal
          integer nlines, nkep, header_length, counter 
          logical initial_growth
    !                                                                       
          character*1024 filin,filout,eos_table
          character*1 ajunk       
    !                                                                                                                           
          mcut=40. 
          max=0                                                    
    !                                                                       
    !--read options                                                         
    !                                                                       
          open(521,file='setup_prep') 
          read(521,*) 
          read(521,522) filin 
          read(521,*) 
          read(521,*) 
          read(521,522) filout 
          read(521,*) 
          read(521,*) 
          read(521,*) grid_goal
          read(521,*) 
          read(521,*)                     
          read(521,*) conv_grid_goal
          read(521,*) 
          read(521,*) 
          read(521,*) enclmass_conv_cutoff
          read(521,*)
          read(521,*)
          read(521,*) deltam(1) 
          read(521,*)
          read(521,*) 
          read(521,*) deltam_growth
          read(521,*)
          read(521,*) 
          read(521,*) maxrad 
          read(521,*) 
          read(521,*) 
          read(521,*) ieos      
          read(521,*) 
          read(521,*)       
          read(521,522) eos_table      
          maxrad = maxrad / udist 
      522 format(A) 
    !
    !--check whether sukhbold's or woosley's progenitors are used
    !
          if (verify('woosley',filin).eq.0) then
              header_length=22
          elseif (verify('sukhbold',filin).eq.0) then
              header_length=2
          else
              print*, "ERROR: unknown dataset is used. Input should have either 'woosley' or 'sukhbold' in its path."
              call exit
          endif
    !                                                                       
    !--get number of entries                                                
    !                                                                       
          nlines = 0 
          OPEN (1, file = trim(filin)) 
          DO 
              READ (1,*, END=42) 
              nlines = nlines + 1 
          END DO 
       42 CLOSE (1) 
    !                                                                       
          nkep = nlines-header_length 
          allocate(vel(nkep)) 
          allocate(rad(nkep)) 
          allocate(dens(nkep)) 
          allocate(t9(nkep)) 
          allocate(yel(nkep)) 
          allocate(ab(nkep)) 
          allocate(omega(nkep)) 
          allocate(press(nkep)) 
    !                                                                       
    !--read-in the data                                                     
    !                                                                       
          open(11,file=trim(filin)) 
          do i=1,header_length 
             read(11,*)ajunk 
          end do 
          do 10 i=1,nkep 
              if (verify('woosley',filin).eq.0) then
                read(11,*) izone, djunk, djunk, drad, dvel, ddens, dtemp,      &
                    djunk,djunk,djunk,domega,dabar,dyel,(yccin(i,j),j=1,19)
              elseif (verify('sukhbold',filin).eq.0) then
                read(11,*) izone, djunk, drad, dvel, ddens, dtemp,             &
                    djunk,djunk,djunk,domega,dabar,dyel,(yccin(i,j),j=1,19)
              endif                                               
             rnorm=0. 
             do j=1,19 
                rnorm=rnorm+yccin(i,j) 
             end do 
             do j=1,19 
                yccin(i,j)=yccin(i,j)/rnorm 
             end do 
    !--rescale to correct units                                             
             rad(i)=drad/udist 
             dens(i)=ddens/udens 
             t9(i)=dtemp/utemp 
             yel(i)=dyel 
             ab(i)=dabar 
             vel(i)=dvel/uvel 
             omega(i)=0. 
       10 continue 
          close (11) 
          dmtot=0. 
          rold=0. 
    
          do i=1,nkep 
             dmtot=dmtot+4.d0*3.14159258/3.d0*dens(i)*                      &
                  (rad(i)**3-rold**3)                                                                                          
             rold=rad(i)                 
          end do 
    
          print *, 'dmtot, rold:      ', dmtot,rold 
          print *, 'initial cell mass:', deltam(1) 
    !                                                                       
    !--nucdata is in nse5.f                                                 
    !                                                                       
          call nucdata 
    !                                                                       
    !--loadmx is in sleos.f                                                 
    !  
          if (ieos.eq.4) then                                                                            
            call loadmx () 
          endif        
    !                   
    !--load eos table
    !      
          if (ieos.eq.5) then 
            call readtable(trim(eos_table))
          endif   
                    
    !
    !--find the optimal deltam for the given convective region resolution      
    !
          call search_deltam(nkep,yccin,vel,rad,dens,              &
                             t9,yel,ab,omega,press,                &
                             maxrad,deltam_growth,ieos,            &
                             conv_grid,conv_grid_goal,             &
                             enclmass_conv_cutoff)
    !
    !--interpolate and setup the grid
    !--growth rate past 200km will keep adjusting 
    !--via a binary-search-like algorithm until 
    !--total grid size is within 1% of the idim
    !     
          call search_deltam_growth(nkep,yccin,vel,rad,dens,       &
                             t9,yel,ab,omega,press,                &
                             maxrad,deltam_growth,ieos,            &
                             conv_grid,enclmass_conv_cutoff,grid_goal)
    !
    !--if you want to avoid binary search and do the original grid setup
    !--comment the above 2 subroutine calls, and uncomment setup_grid() below
    !                             
        !   call setup_grid(nkep,yccin,vel,rad,dens,                &
        !                   t9,yel,ab,omega,press,                  &
        !                   maxrad,deltam_growth,ieos,              &
        !                   conv_grid,enclmass_conv_cutoff)          
                                  
      102 format(20(1pe9.2)) 
      103 format(3(1pe16.8)) 
      115 format(I5,1pe13.4) 
          open (29,file=trim(filout),form='unformatted') 
          call wdump 
          print *,'--------------- Summary ---------------'
          print *,'  Input progenitor:          ', trim(filin)
          print*, '  Final initial cell mass:', deltam(1)          
          print *,'  Final growth rate:      ', deltam_growth 
          print *,'  Convective Grid Size:', conv_grid         
          print *,'  Grid Size:           ', ncell 
          print *,'  Output file name:            ', trim(filout) 
          print *,'---------------------------------------' 
          stop 
          END                                                                                     
!      
  
          subroutine setup_grid(nkep,yccin,vel,rad,dens,                   &
                                t9,yel,ab,omega,press,                     &
                                maxrad,deltam_growth,ieos,                 &
                                conv_grid,enclmass_conv_cutoff)                         
!****************************************************************              
!                                                               *              
! Woosely's code has nkep cells.  I want to make my code        *       
! whatever number of cells I like.  I use an average of         *      
! adjacent cells to calculate the important values of my code.  *       
!                                                               *              
!****************************************************************                                  
          implicit double precision (a-h, o-z) 
          !                                                                       
          parameter(utemp=1e9) 
          parameter(udens=2e6) 
          parameter(uvel=1e8) 
          parameter(udist=1e9) 
          parameter(utime=1e1) 
          parameter(uergg=1e16) 
          parameter(pi43=3.14159*4.0/3.0) 
          !                                                                       
          parameter(idim=4000) 
          !                                                                       
          common /celle/ x(0:idim),v(0:idim) 
          common /cellc/ u(idim),rho(idim),ye(idim),q(idim),dq(idim) 
          common /numb/ ncell 
          common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),        &
          &               unue(idim),unueb(idim),unux(idim)                  
          common /eosq / pr(idim), vsound(idim), u2(idim), vsmax 
          common /carac/ deltam(idim), abar(idim) 
          common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
          common /freez/ ufreez(idim) 
          common /tempe/ temp(idim) 
          common /cent/ dj(idim) 
          common /abun/ ycc 
  
          dimension vel(nkep),rad(nkep),dens(nkep),t9(nkep),        &
          yel(nkep),ab(nkep),omega(nkep),press(nkep) 
          dimension yccin(idim,20)
  
          real ycc(idim,20) 
          double precision rhocgs, tkelv, yej, abarj, ucgs, pcgs, xpj, xnj 
          double precision enclmass(0:idim)                       
          double precision maxrad, deltam_growth         
  
          logical beyond_focus
          integer nkep,ieos,conv_grid,one_counter,pns_counter
               
          conv_grid = 0
          one_counter = 0
          beyond_focus = .false.         
          pns_counter = 0                                                        
          x(0) = 0.d0 
          v(0) = 0.d0 
          enclmass(0) = 0.d0 
  
          do i=1,idim 
            do k=1,nkep 
              if (rad(k).gt.x(i-1)) goto 20 
            end do 
     20    continue 
  
           if (k.eq.1) then !-- for k=1, we can not use cell k-1                         
              ind      = k
              fac      = x(i-1)/rad(1) 
              fac1     = 1-fac
              fac2     = fac
              facomei1 = fac1
              facomei2 = fac2
           else !--for all other k      
              ind      = k-1
              fac      =-(x(i-1)-rad(k))/(rad(k)-rad(k-1)) 
              fac1     = fac
              fac2     = 1-fac
              facomei1 = fac2
              facomei2 = fac1
           endif
  
          rho(i)  = fac1*0.5*dens(ind) +                              &
                        0.5*dens(k)+fac2*0.5*dens(k+1)
          v(i)    = fac1*0.5*vel(ind) +                               &
                        0.5*vel(k)+fac2*0.5*vel(k+1)
          temp(i) = fac1*0.5*t9(ind) +                                &
                        0.5*t9(k)+fac2*0.5*t9(k+1)
          ye(i)   = fac1*0.5*yel(ind) +                               &
                        0.5*yel(k)+fac2*0.5*yel(k+1)
          abar(i) = fac1*0.5*ab(ind) +                                &
                        0.5*ab(k)+fac2*0.5*ab(k+1)
          x(i)    = (x(i-1)**3+deltam(i)/(pi43*rho(i)))**0.3333333333 
          
          do j=1,19 
              ycc(i,j)=fac1*0.5*yccin(ind,j) +                        &
      &           0.5*yccin(k,j)+fac2*0.5*yccin(k+1,j)
          end do 
          
          if (k.eq.1) then
              ycc(i,3)=ycc(i,1) 
              do j=6,15 
                ycc(i,j)=ycc(i,j+1) 
              end do 
              ycc(i,16)=ycc(i,18) 
              ycc(i,17)=ycc(i,19)
          endif
                                                                          
          omei=facomei1*0.5*omega(k) +                                &
                 0.5*omega(k)+facomei2*0.5*omega(k+1)                        
          dj(i)=omei*x(i)*x(i) 
          
          enclmass(i) = enclmass(i-1) + deltam(i) 
  !                                                                        
  !--Adjust the grid here, choose a focus region
  !
          write(45,*) i,enclmass(i),dj(i) 
      
          grid_size = int(idim)
          if (enclmass(i).le.enclmass_conv_cutoff) then
              deltam(i+1) = deltam(1)
          else
              if (beyond_focus.eqv..false.) then
                  leftover = grid_size-i             
                  beyond_focus = .true.
                  focus_i = i
              endif

              deltam(i+1) = deltam(focus_i)*(1+deltam_growth)**(i-focus_i)

          end if 
          if (k.eq.1) one_counter = one_counter+1
  !                                                                       
  !--call eos                                                             
  !                                                                       
          rhocgs=dble(rho(i)*udens) 
          tkelv=dble(temp(i)*utemp) 
          yej=dble(ye(i)) 
          abarj=dble(abar(i)) 
          iflag=0 
  
          ! choose either Swesty-Lattimer or SFHo Tables EOS
          if (ieos.eq.4) then
              call eosgen(i,iflag,rhocgs,tkelv,yej,abarj,                 &
                          ucgs,pcgs,xpj,xnj,scgs)                                
          elseif (ieos.eq.5) then
              call eos5(i,rho,u,temp,ye,rhocgs,tkelv,yej,            &
                          abarj,ucgs,pcgs,xpj,xnj,scgs,iflag)       
          endif

          u(i)=ucgs/uergg 
          pr(i)=pcgs/uergg/udens 
          u2(i)=scgs/uergg*utemp 

          if (iflag.eq.1.) then 
              if (abarj.lt.2.5) then 
                  ufr=-3.3d17/uergg 
              elseif (abarj.lt.4.5) then 
                  ufr=-6.1d17/uergg 
              elseif (abarj.lt.20.) then 
                  ufr=-7.7d18/uergg 
              else 
                  ufr=-8.4d18/uergg 
              endif 
          else 
              ufr=0.0 
          endif 

          ufreez(i)=ufr 
          ifleos(i)=iflag 
          xp(i)=xpj 
          xn(i)=xnj 
          
          !if (x(i).le.2.d-2) conv_grid = conv_grid+1
          if (enclmass(i).le.enclmass_conv_cutoff) conv_grid = conv_grid+1

          if (x(i).gt.maxrad) then 
              ncell=i 
              goto 50 
          end if
        enddo 
        ! print*, enclmass(i-10:i)
        ! print*, 'max', maxval(enclmass)
        ! stop
  !                                                                               
        50 continue 

         return
         end        


         subroutine search_deltam(nkep,yccin,vel,rad,dens,              &
                                  t9,yel,ab,omega,press,                &
                                  maxrad,deltam_growth,ieos,            &
                                  conv_grid,conv_grid_goal,             &
                                  enclmass_conv_cutoff)                     
!*********************************************************              
!                                                        *              
!  Find the optimal deltam for the given convective      *
!  region resolution                                     *             
!                                                        *              
!*********************************************************  

         implicit double precision (a-h, o-z) 

         parameter(idim=4000)
         dimension vel(nkep),rad(nkep),dens(nkep),t9(nkep),            &
                   yel(nkep),ab(nkep),omega(nkep),press(nkep) 
         dimension yccin(idim,20)
         
         common /carac/ deltam(idim), abar(idim)
 
         double precision enclmass(0:idim)                       
         double precision maxrad, deltam_growth  

         integer ieos
         integer conv_grid, conv_grid_goal
         integer nkep, counter 
         logical initial_growth

         print *,'-------- Search for deltam(1) ---------'
         initial_growth = .true.
         adjust_deltam = 1.d1         
         counter = 1
         low = 0.0
         
         do while (.true.)         
 
             call setup_grid(nkep,yccin,vel,rad,dens,                   &
                             t9,yel,ab,omega,press,                     &
                             maxrad,deltam_growth,ieos,                 &
                             conv_grid, enclmass_conv_cutoff) 

             if (conv_grid.eq.0) then
                 print*, "ERROR: Initial Cell Mass is to small for the grid! Increase it."
                 call exit
             else
                 write(*,'(I5,A,I5)'), counter, ' Convective size:', conv_grid
             endif
 
             ! exit condition once convection grid is of the right size
             if (conv_grid.eq.conv_grid_goal) then
                 exit
             
             elseif (conv_grid.lt.conv_grid_goal) then            
                 if (.not.initial_growth) then
                     low = adjust_deltam
                     deltam(1) = deltam(1)*adjust_deltam           
                     adjust_deltam = adjust_deltam + (high-adjust_deltam)*0.5                
                 endif
                 deltam(1) = deltam(1)/adjust_deltam
 
             elseif (conv_grid.gt.conv_grid_goal) then
                 if (initial_growth) initial_growth = .false.
                 high = adjust_deltam
                 deltam(1) = deltam(1)*adjust_deltam
                 adjust_deltam = adjust_deltam - (adjust_deltam-low)*0.5
                 deltam(1) = deltam(1)/adjust_deltam
             endif

             counter = counter + 1

             if (counter .eq. 999) then
               print*, 'ERROR: Could not converge on requested convection grid size'
               call exit
             endif
         enddo
         print *,'---------------------------------------'
         write(*,'(A,I3,A)'), ' Cell mass in convective region converged in', counter-1, ' iterations'
         print*,'Final deltam(1): ', deltam(1)
         print*, ''
      end

      subroutine search_deltam_growth(nkep,yccin,vel,rad,dens,       &
                               t9,yel,ab,omega,press,                &
                               maxrad,deltam_growth,ieos,            &
                               conv_grid,enclmass_conv_cutoff,grid_goal)                     
!*********************************************************              
!                                                        *              
!  Growth rate past 200km will keep adjusting            * 
!  via a binary-search-like algorithm until              *
!  total grid size is within 1% of the idim              *          
!                                                        *              
!*********************************************************  

      implicit double precision (a-h, o-z) 
  
      parameter(idim=4000)
      dimension vel(nkep),rad(nkep),dens(nkep),t9(nkep),            &
      yel(nkep),ab(nkep),omega(nkep),press(nkep) 
      dimension yccin(idim,20)
  
      common /numb/ ncell
      common /carac/ deltam(idim), abar(idim)
  
      double precision enclmass(0:idim)                       
      double precision maxrad, deltam_growth  
  
      integer ieos
      integer conv_grid,grid_goal
      integer nkep, counter 
      logical initial_growth

      print *,'---- Search for deltam growth rate ----'
      initial_growth = .true.
      adjust_growth = 1.d1
      counter = 1
      low = 0.0 
      high = -1

      do while (.true.) 
          ncell = 0   

          call setup_grid(nkep,yccin,vel,rad,dens,                   &
                          t9,yel,ab,omega,press,                     &
                          maxrad,deltam_growth,ieos,                 &
                          conv_grid,enclmass_conv_cutoff) 
                        
          ! exit condition once precision of under 1% is reached
          if (abs(grid_goal-ncell).lt.0.01*ncell) then
            exit

          ! increase growth rate if the grid is above grid_goal
          elseif (ncell.eq.0.or.ncell.gt.grid_goal) then
              
              write(*,'(I5,A,I5,A)'), counter,' Grid size: above grid_goal! (', grid_goal,')'

              if (.not.initial_growth) then
                  low = adjust_growth
                  deltam_growth = deltam_growth/adjust_growth
                  adjust_growth = adjust_growth + (high-adjust_growth)*0.5
              endif

              deltam_growth = deltam_growth*adjust_growth

          ! decrease growth rate if the grid is below grid_goal   
          else
              write(*,'(I5,A,I5)'), counter, ' Grid size:', ncell
              
              if (initial_growth) initial_growth = .false.
              high = adjust_growth
              deltam_growth = deltam_growth/adjust_growth            
              adjust_growth = adjust_growth - (adjust_growth-low)*0.5
              deltam_growth = deltam_growth*adjust_growth
          endif          

          counter = counter + 1
          
          if (counter .eq. 999) then
              print*, 'ERROR: Could not converge on deltam_growth.'
              print*, '       Probably not enough mass enclosed in the maximum radius region.'
              call exit
          endif
      enddo
      print *,'---------------------------------------'
      write(*,'(A,I3,A)'), ' Growth rate converged in', counter, ' iterations'
      print*, 'Final deltam_growth:', deltam_growth   
      print*, ''   

      end
        
      subroutine wdump 
!************************************************************           
!                                                           *           
!  this routine writes a dump on disk                       *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h, o-z) 
!                                                                       
      logical from_dump 
      integer idump 
      parameter (idim=4000) 
      common /celle/ x(0:idim),v(0:idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim),dq(idim) 
      common /numb/ ncell 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)                  
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /eosq / pr(idim), vsound(idim), u2(idim), vsmax 
      common /freez/ ufreez(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /cgas / gamma 
      common /times / t, dt 
      common /ener2/ tkin, tterm 
      logical te(idim), teb(idim), tx(idim) 
      common /cent/ dj(idim) 
      real ycc(idim,20) 
      common /abun/ ycc 
      common /timei/ steps(idim) 
      common /rshock/ shock_ind, shock_x 
      common /pns/ pns_ind, pns_x 
      common /dump/ from_dump 
      double precision rlumnue,rlumnueb,rlumnux,bounce_time 
      double precision pr_turb(idim) 
!            
      steps = 0 
      shock_ind = 0 
      shock_x = 0 
      pns_ind = 0 
      pns_x = 0 
      idump=0 
      bounce_time=0 
      from_dump = .false. 
!--initialize neutrino fluxes                                           
      rlumnue = 0 
      rlumnueb = 0 
      rlumnux = 0 
      pr_turb(:) = 0 
!                      
      nc = ncell 
      do i=1,nc 
         if (u(i).eq.0) then
            print*, "ERROR: u(i) is 0 at i=",i
            stop 
         endif
         ynue(i)=0. 
         ynueb(i)=0. 
         ynux(i)=0. 
         unue(i)=0. 
         unueb(i)=0. 
         unux(i)=0. 
         te(i)=.false. 
         teb(i)=.false. 
         tx(i)=.false. 
         q(i)=0. 
         dq(i)=0. 
      enddo 
      rb=9e-3 
      gc=0.0001 
      ftrap=1. 
      fe=ftrap 
      fb=ftrap 
      fx=ftrap 
!                                                                       
!--write                                                                
!                     
      nqn=17 
      write(29,iostat=io,err=10) idump,nc,t,gc,rb,fe,fb,fx,             &
     &     pns_ind,pns_x,shock_ind,shock_x,                             &
     &     bounce_time,from_dump,rlumnue,rlumnueb,rlumnux,              &
     &     (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),    &
     &     (u(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),           &
     &     (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),             &
     &     (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),            &
     &     (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),         &
     &     (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),         &
     &     (ufreez(i),i=1,nc),(pr(i),i=1,nc),(u2(i),i=1,nc),            &
     &     (dj(i),i=1,nc),                                              &
     &     (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),               &
     &     (steps(i),i=1,nc),((ycc(i,j),j=1,nqn),i=1,nc),               &
     &     (vsound(i),i=1,nc),(pr_turb(i),i=1,nc)                       
!                                                     
     !print*, 'abar max & min: ', maxval(abar), minval(abar)
    !  print*, temp
    !  print*, 'temp:', temp(nc-10:nc)
    !  print*, 'abar:', abar(nc-10:nc)
      do i=1,ncell 
         write (43,103) (ycc(i,j),j=1,19) 
      end do 
!      print *, x(0),x(1)                                               
  102 format(I4,4(1pe14.4),I3) 
  103 format(19(1pe12.4)) 
      return 
!                                                                       
!--an error as occured while writting                                   
!                                                                       
   10 print *,'an error has occured while writing' 
!                                                                       
      return 
      END                                           
      
      subroutine eos5(i,rho,u,temp,ye_table,rhocgs,tkelv,yej,            &
                    abarj,ucgs,pcgs,xpj,xnj,scgs,iflag)

      use eosmodule
      
      implicit double precision (a-h,o-z)

      parameter(idim=4000)
      parameter(utemp=1e9) 
      parameter(udens=2e6) 
      parameter(uvel=1e8) 
      parameter(udist=1e9) 
      parameter(utime=1e1) 
      parameter(uergg=1e16) 
      parameter(pi43=3.14159*4.0/3.0)  
!     
      dimension u(idim),rho(idim),ye_table(idim),temp(idim)

      data ggcgs/6.67e-8/, avo/6.02e23/ 
      data aradcgs /7.565e-15/, boltzk/1.381e-16/ 
      data hbar/1.055e-27/, ccgs/3e10/ 
      data emssmev/0.511e0/, boltzmev/8.617e-11/ 
      data ergmev /6.2422e5/, sigma1/9d-44/, sigma2/5.6d-45/ 
      data c2cgs /6.15e-4/, c3cgs /5.04e-10/, fermi/1d-13/

      real*8 xrho,xye,xtemp,xtemp2
      real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
      real*8 xdpderho,xdpdrhoe
      real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat_table
      real*8 xxa,xxh,xxn,xxp              
      integer keytemp,keyerr 
      
      keytemp = 1
      keyerr  = 0

      xrho  = rho(i)*udens
      xenr  = u(i)*uergg
      xtemp = temp(i)*utemp*boltzmev
      xye   = ye_table(i)

      if (xye  .ge.0.6  ) then
          print*, "WARNING: ye > yemax, setting ye to maxye at", i
          xye  = 0.599
      endif
      if (xrho .lt.166.5) then
          print*, "WARNING: rho < rhomin, setting rho to rhomin at", i
          xrho = 166.5 
      endif
      if (xtemp.lt.0.011) then
          print*, "WARNING: temp < tempmin, setting temp to tempmin at", i
          xtemp= 0.011                                
      endif
      !                                                                       
      !--SFHo EOS tables                                               
      !                     
      call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,               &
            xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,    &
            xmuhat_table,keytemp,keyerr,precision)  

      rhocgs = xrho
      tkelv  = xtemp/boltzmev
      yej    = xye
      abarj  = xabar
      ucgs   = xenr
      pcgs   = xprs
      xpj    = xxp
      xnj    = xxn
      scgs   = xent  
      iflag  = 5 
      end

      subroutine eosgen(izone,iflag,rhocgs,tkelv,ye,abar,               &
     &                  ucgs,pcgs,xp,xn,scgs)                           
                                                                        
      implicit double precision (a-h,o-z) 
                                                                        
      parameter(t9nse=7d9) 
      parameter(rhoswe=1.d11) 
                                                                        
      double precision inpvar(4) 
      common /pvar/ inpvar 
      if (iflag.eq.0) then 
         if (tkelv.lt.t9nse.and.rhocgs.lt.rhoswe) then 
             iflag=1 
         else 
            if (rhocgs.lt.rhoswe) then 
               iflag=2 
            else 
               iflag=3 
               print *, 'sleos' 
            endif 
         endif 
      endif 
      if (iflag.eq.1) then 
         zbar=ye*abar 
         call coulomb(rhocgs,zbar,ye,ucoul,pcoul) 
         tno=tkelv/1d9 
         rhono=rhocgs/1d7 
         call nados(tno,rhono,zbar,abar,pel,eel,sel,                    &
     &              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)            
         ucgs=etot*1d17+ucoul 
         pcgs=ptot*1d24+pcoul 
         scgs=stot*1d8 
         xp=0.0 
         xn=0.0 
      elseif (iflag.eq.2) then 
         t9=tkelv/1d9 
         call nsestart(t9,rhocgs,ye,xp,xn) 
         call nsetemp(ipart,t9,rhocgs,ye,t9,xp,xn,                      &
     &                zbar,abar,ubind,dubind)                           
         zbar=ye*abar 
         call coulomb(rhocgs,zbar,ye,ucoul,pcoul) 
         tno=tkelv/1d9 
         rhono=rhocgs/1d7 
         call nados(tno,rhono,zbar,abar,pel,eel,sel,                    &
     &              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)            
         ucgs=etot*1d17+ucoul+ubind 
         pcgs=ptot*1d24+pcoul 
         scgs=stot*1d8 
      elseif (iflag.eq.3) then 
         tswe=tkelv/1.16d10 
         inpvar(1)=tswe 
         inpvar(2)=0.155d0 
         inpvar(3)=-15.0d0 
         inpvar(4)=-10.d0 
         brydns=rhocgs*6.02d-16 
         pprev=ye*brydns 
!        print *,'calling slwrap: ye',ye,brydns,tkelv,tswe              
         call slwrap(ipart,inpvar,ye,brydns,pprev,                      &
     &               psl,usl,dusl,gamsl,eta,xp,xn,ssl)                  
         abar=1.d0 
         ucgs=usl*9.644d18 
         pcgs=psl*1.602d33 
         scgs=ssl*8.25d7 
      endif 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine coulomb(rhoi,zbar,ye,ucoul,pcoul) 
!***********************************************************            
!                                                                       
!  compute Coulomb corrections as given in Shapiro                      
!  and Teukolsky. p. 31 (2.4.9) and (2.4.11)                            
!  in cgs:                                                              
!        ucoul=-1.45079*e**2*avo**4/3*ye**4/3*rho**1/3*Z**2/3           
!             =-1.70e13 Ye**4/3 * rho**1/3 * Z**2/3                     
!  code units: mulitply by udens**1/3 / uergg                           
!                                                                       
!        pcoul=-0.4836*e**2*avo**4/3*Ye**4/3*rho**4/3*Z**2/3            
!             =-5.67e12 Ye**4/3 * rho**4/3 * Z**2/3                     
!  code units: mulitply by udens**4/3 / uergcc                          
!                                                                       
!***********************************************************            
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
!     parameter(ufac=-0.214d0)                                          
!     parameter(pfac=-0.0714d0)                                         
      parameter(ufac=-1.70d13) 
      parameter(pfac=-5.67d12) 
!                                                                       
      rho13=rhoi**0.333333333333d0 
      rho43=rho13*rhoi 
      ye2=ye*ye 
      y43z23=(ye2*zbar)**0.66666666666d0 
      ucoul=ufac*rho13*y43z23 
      pcoul=pfac*rho43*y43z23 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine slwrap(ipart,inpvar,yesl,brydns,pprev,                 &
     &                  psl,usl,dusl,gamsl,etasl,ypsl,ynsl,ssl)         
!******************************************************************     
!                                                                       
!  This is the wrapper routine for the swesty-lattimer eos.             
!                                                                       
!******************************************************************     
!                                                                       
!-- 0.46*(mn-mp-me) + bindfe56 MeV/nucleon                              
      double precision ushift 
!     parameter(ushift=0.46d0*0.783d0+8.7904d0)                         
      parameter(ushift=0.0d0) 
!                                                                       
      include 'eos_m4a.inc' 
      include 'el_eos.inc' 
!                                                                       
!--these variables are needed but not declared in the include files     
      integer sf 
      double precision told, pprev, ssl 
      double precision yesl, psl, usl, dusl, gamsl, etasl, ypsl, ynsl 
!                                                                       
!     print *,'slwrap: calling inveos, ipart=',                         
!    1        ipart,inpvar(1),brydns,yesl                               
      ye=yesl 
      call inveos(inpvar,told,ye,brydns,1,eosflg,0,sf,                  &
     &            xprev,pprev)                                          
!                                                                       
      if (sf.ne.1) print *,'inveos fails for particle',ipart 
      psl=ptot 
      usl=utot+ushift 
      dusl=dudt 
!                                                                       
      gamsl=gam_s 
      etasl=musube/inpvar(1) 
!                                                                       
!-- free (exterior) nucleon fractions                                   
      ypsl=xprot 
      ynsl=xnut 
      ssl=stot 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nsetemp(ipart,t9old,rho,ye,t9,yp,yn,                   &
     &                   zbar,abar,ubind,dubind)                        
!*************************************************************          
!                                                                       
! this subroutine figures out the NSE eq. assuming that yp              
! and yn were previously know at the SAME density and ye,               
! but different temperature                                             
!                                                                       
!**************************************************************         
!write(*,*)                                                             
!                                                                       
      implicit double precision(a-h,o-z) 
      parameter (tolnse=1d-5,kmax=10) 
!                                                                       
      common /testnse/ testk,testzy,testay,testyp,testyn 
!                                                                       
!--finding zero point of binding energy                                 
!                                                                       
      ider=2 
      call nsesolv(ider,t9old,rho,ye,yp,yn,kit,kmax,ubind0,             &
     &             zbar,abar)                                           
!                                                                       
      If(kit.ge.kmax) Then 
          write(*,*) 'NSE mis-stored entering nsetemp' 
          write(*,*) 'T9, rho, ye',t9old,rho,ye 
          write(*,*) 'inconsistent with yp, yn',yp,yn 
      Endif 
      ypold=yp 
      ynold=yn 
      t9min=min(t9,t9old) 
      ider=0 
!                                                                       
!--pick initial delt9                                                   
!                                                                       
      delt9=0.5d0 
      if(t9min.lt.12.d0)then 
         delt9=0.05d0 
      elseif(t9min.gt.20.d0) then 
         delt9=2.d0 
      end if 
      delt9=dsign(delt9,t9-t9old) 
      ypold=yp 
      ynold=yn 
      t9last=t9old 
!                                                                       
!--Begin temp iteration                                                 
!                                                                       
      do i=1,2000 
          delt9=dsign(min(dabs(t9-t9last),dabs(delt9)),delt9) 
          t9tmp=t9last+delt9 
          call nsesolv(ider,t9tmp,rho,ye,yp,yn,kit,kmax,ubind,          &
     &                   zbar,abar)                                     
          If(dabs((t9tmp-t9)/t9).lt.tolnse.and.kit.lt.kmax) goto 70 
          If (kit.ge.kmax) Then 
              delt9=.5d0*delt9 
              yp=ypold 
              yn=ynold 
          Elseif(kit.lt.4) Then 
              delt9=2.d0*delt9 
              t9last=t9tmp 
              ypold=yp 
              ynold=yn 
          Else 
              t9last=t9tmp 
              ypold=yp 
              ynold=yn 
          Endif 
      enddo 
!                                                                       
!--did not converge, print out error                                    
!                                                                       
      write(*,*)'nsetemp(1) did not converge!',ipart 
      write(*,*)t9old,t9tmp,t9 
      write(*,*)i,delt9,rho 
      write(*,*)ye,yp,yn 
      write(*,*)kit,kmax,ubind,dubind 
      write(*,*)abar,zbar 
      write(*,*)'before nsestart',yp,yn 
      call nsestart(t9,rho,ye,yp,yn) 
      write(*,*)'after nsestart',yp,yn 
   70 Continue 
!                                                                       
!--solve for actual t9, calculate binding energy, average A, Z          
!                                                                       
      ider=2 
      call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,                 &
     &                   zbar,abar)                                     
      If(kit.ge.kmax) Then 
          write(*,*) 'NSEtemp failed for final T9, particle',ipart 
          write(*,*) kit,t9,t9tmp 
          STOP 
      Endif 
!                                                                       
!--Overstep in T9 to calculate dUb/dT9                                  
!                                                                       
      If(t9.lt.12.d0) Then 
          delt9=1d-4 
      Elseif(t9.gt.20.d0) Then 
          delt9=1.d0 
      Else 
          delt9=1d-2 
      Endif 
   80 t9d=t9+delt9 
      ypd=yp 
      ynd=yn 
      call nsesolv(ider,t9d,rho,ye,ypd,ynd,kit,kmax,ubindd,             &
     &             zbar,abar)                                           
      If(kit.ge.kmax) Then 
          write(*,*) 'dUb/dT step failure in nsetemp, particle ',ipart 
          write(*,*) kit,t9d,t9 
          delt9=.5d0*delt9 
          goto 80 
      Endif 
      dubind=(ubindd-ubind)/(t9d-t9) 
!                                                                       
!                                                                       
!--step back in T9, in order to calculate dUbind/dT9                    
!                                                                       
!     call nsesolv(ider,t9last,rho,ye,ypold,ynold,kit,kmax,ubindlast,   
!    &                   zbar,abar)                                     
!     If(kit.ge.kmax) Then                                              
!         write(*,*) 'NSEtemp failed for deriv T9last: particle',ipart  
!         write(*,*) kit,t9,t9tmp                                       
!         STOP                                                          
!     Endif                                                             
!     dubind=(ubindlast-ubind)/(t9last-t9)                              
      Return 
      END                                           
!                                                                       
                                                                        
      subroutine nsestart(t9fnl,rhofnl,yefnl,yp,yn) 
      parameter (nsetol=1d-4,kmax=7) 
      implicit double precision (a-h,o-z) 
      T9start=100.d0 
      rho=rhofnl 
      ye=yefnl 
      yp=ye 
      yn=1.d0-ye 
      ider=0 
      T9=T9start 
      DelT9=dsign(20.d0,(t9fnl-t9start)) 
  100 testt=t9fnl-t9 
!     write(*,*) 'testt',testt                                          
      If(dabs(testt).gt.nsetol) Then 
          call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,             &
     &                   zbar,abar)                                     
!         write(*,999) t9fnl,t9,rho,ye,yp,yn,kit                        
          If(kit.ge.kmax) Then 
              t9=t9-delt9 
              delt9=.5d0*delt9 
          Else 
              t9last=t9 
              If(kit.le.3) Then 
                 delt9=2.d0*delt9 
              Endif 
          Endif 
          If(t9.gt.35.d0) Then 
              dtmax=10.d0 
          Elseif(t9.gt.25.d0) Then 
              dtmax=3.d0 
          Elseif(t9.gt.10.d0) Then 
              dtmax=1.d0 
          Elseif(t9.gt.4.5d0) Then 
              dtmax=.1d0 
          Else 
              dtmax=.01d0 
          Endif 
          delT9=dsign(dmin1(dabs(testt),dabs(delt9),dabs(dtmax)),testt) 
          T9=T9+delT9 
          GOTO 100 
      Else 
          t9=t9fnl 
          ider=2 
          kmtmp=3*kmax 
          call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmtmp,ubind,            &
     &                   zbar,abar)                                     
!         write(*,999) t9fnl,t9,rho,ye,yp,yn,kit                        
          If(kit.ge.kmtmp) Then 
             print*,'broken code' 
             stop 
          Endif 
      Endif 
      Return 
  999 Format(6(1pd11.4),1x,I3) 
                                                                        
      END                                           
                                                                        
                                                                        
      subroutine nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,           &
     &                   zbar,abar)                                     
      parameter (nyab= 16, nyabmx= 16) 
      implicit double precision (a-h,o-z) 
      dimension y(nyabmx),dlg(24,nyabmx),x(nyabmx) 
      dimension za(nyabmx),zp(nyabmx),zn(nyabmx),t9i(24) 
      dimension h(37),dhdp(37),dhdn(37),dhdt9(37),intz(nyabmx) 
      character*5 inam(nyabmx) 
      common /zanp/ za,zp,zn,pnp,bi(nyabmx),angm(nyabmx) 
      common /alph/ dlgalp(nyabmx) 
      common /name/ inam 
      common /partf/ dlg,t9i 
      common /testnse/ testk,testzy,testay,testyp,testyn 
      iscrn=2 
      If(ider.eq.2) Then 
          tol=1d-6 
      Else 
          tol=1d-4 
      Endif 
      If(rho.gt.1d11) Then 
          T9npa=50.d0 
      Elseif(rho.gt.1d10) Then 
          T9npa=25.d0 
      Elseif(rho.gt.1d9) Then 
          T9npa=16.d0 
      Elseif(rho.gt.1d8) Then 
          T9npa=12.d0 
      Elseif(rho.gt.1d7) Then 
          T9npa=10.d0 
      Else 
          T9npa=8.d0 
      Endif 
      If(t9.gt.t9npa) Then 
          nmax=3 
      Else 
          nmax=nyab 
      Endif 
      yp0=yp 
      yn0=yn 
!     write(*,*) 'NSEsolv',t9,rho,ye                                    
!--Raph's routine to find abar, zbar, and nuclear binding energy        
!  ---------------------------------------------------------------------
!  This portion is designed to calculate the Nuclear Statistical equilib
!  distribution of material for a given temperature and density.  Analyt
!  expressions exist which give solutions for Y(A,Z) as a product of sep
!  functions of Yp, Yn, and alpha, where alpha contains all of the tempe
!  density, and nuclear parameter dependence.  Yp and Yn are the free pr
!  and neutron abundances.  Implicit solution of equations for Ye and ma
!  conservation, using the analytic expressions to convert Y(A,Z) to Yp 
!  yields values for Yp and Yn and hence the entire distribution.       
!  We will employ a two dimensional Newton method, which, though it conv
!  extremely rapidly if it is near the solution, is prone to severe     
!  convergence problems.  For this reason the subroutine uses the values
!  Yp and Yn from the previous step.  The subroutine also reports if it 
!  converges.  In the event that it fails to converge, the correct solut
!  requires a more detailed iteration in T, rho and Ye.                 
!  ---------------------------------------------------------------------
      dlgyp=dlog(yp) 
      dlgyn=dlog(yn) 
      Do 100 i=1,nmax 
        If(zp(i).gt.1.d0) Then 
            intz(i)=int(zp(i)) 
        Else 
            intz(i)=1 
        Endif 
  100 Continue 
!--Calculate the thermodynamic dependence coefficients                  
      call alphcalc(nmax,t9,rho) 
!     call nsenorm(y,za,nmax,nmax,zp,zn,ye)                             
!--Calculate the screening factors                                      
      call nsescreen(T9,rho,ye,h,dhdp,dhdn,dhdT9,iscrn) 
      k=1 
      testk=1.d0+tol 
!  ---------------------------------------------------------------------
!  For each iterative T and rho it is necessary to first calculate alph(
!  which contains the dependence on temperature, density and nuclear par
!  In the case of screening (iscrn=2) the alph(n) are then modified by t
!  appropriate screening factor, the product of the screening factors fo
!  the series of (p,gamma) reactions necessary to build up the required 
!  Once alph(n) has been calculated, the next step is to calculate f, th
!  implicit formula for Ye, g, the implicit formula for mass conservatio
!  dfp, dfn, dgn, and dgp, the derivatives of f and g with respect to yp
!  These 6 values, calculated using logaritmic interim steps to avoid ov
!  and the analytic solution of the matrix eqn, yield delyp and delyn. T
!  new values of yp,yn are yp+delyp, yn+delyn. The implicit solution con
!  until either kmax iterations are run, or the change in solution is le
!  than some minimum tolerance, tol.                                    
!  ---------------------------------------------------------------------
   20 If(k.lt.kmax) THEN 
          IF(testk.gt.tol) THEN 
            f=-ye 
            g=-1.d0 
            dfp=0.d0 
            dfn=0.d0 
            dgp=0.d0 
            dgn=0.d0 
!           dfp2=0.d0                                                   
!           dfn2=0.d0                                                   
!           dgp2=0.d0                                                   
!           dgn2=0.d0                                                   
            dlgyp=dlog(yp) 
            dlgyn=dlog(yn) 
            ypm1=1.d0/yp 
            ynm1=1.d0/yn 
            DO 30 n=1,nmax 
              iz=intz(n) 
              zpt=zp(n) 
              znt=zn(n) 
              zat=za(n) 
              yt=dexp(dlgalp(n)+zpt*dlgyp+znt*dlgyn+h(iz)) 
              If(yt.gt.1.d+100) Then 
!                 write(*,*) 'y explosion'                              
                  k=kmax 
                  Goto 20 
              Endif 
              dyp=yt*ypm1 
              dyn=yt*ynm1 
              f=f+zpt*yt 
              g=g+zat*yt 
              zpdyp=zpt*dyp 
              zndyn=znt*dyn 
              dfp=dfp+zpt*zpdyp 
!             dfp2=dfp2+zpt*dhdp(iz)*yt                                 
              dfn=dfn+zpt*zndyn 
!             dfn2=dfn2+zpt*dhdn(iz)*yt                                 
              dgp=dgp+zat*zpdyp 
!             dgp2=dgp2+zat*dhdp(iz)*yt                                 
              dgn=dgn+zat*zndyn 
!             dgn2=dgn2+zat*dhdn(iz)*yt                                 
   30       CONTINUE 
!           write(*,*) 'K',k,f,g                                        
!           dfn=dfn+dfn2                                                
!           dfp=dfp+dfp2                                                
!           dgn=dgn+dgn2                                                
!           dgp=dgp+dgp2                                                
            det=dfp*dgn-dfn*dgp 
            IF(det.ne.0.0d0)THEN 
                ddet=1.d0/det 
                delyp=(dgn*f-dfn*g)*ddet 
                delyn=(dfp*g-dgp*f)*ddet 
            ELSE 
                delyp=0.d0 
                delyn=0.d0 
                k=kmax 
            ENDIF 
            ypt=yp-delyp 
            ynt=yn-delyn 
            If(ypt.gt.0.d0.and.ynt.gt.0.d0) Then 
                yp=ypt 
                yn=ynt 
            Else 
                k=kmax 
            Endif 
!           write(10,4000)k,yp,yn                                       
!           write(10,4300)f,g                                           
!           write(15,4200)dfp,dfn,dgp,dgn                               
!           write(15,4250)dfp2,dfn2,dgp2,dgn2                           
!           write(10,4100)delyp,delyn,det                               
            k=k+1 
            testzy=f 
            testay=g 
            testyp=delyp/yp 
            testyn=delyn/yn 
            testk=sqrt(testyp*testyp+testyn*testyn+f*f+g*g) 
            GOTO 20 
          ENDIF 
!--NSE converges                                                        
          kit=k 
!         write(10,4000)k,yp,yn                                         
!         write(10,3100)T9r,rhor,ye                                     
      ELSE 
!--NSE fails to converge in Kmax steps                                  
!         write(10,3200)kmax,T9r,rhor                                   
          yp=yp0                                                        
          yn=yn0                                                        
          kit=k                                                         
      ENDIF                                                             
!  ---------------------------------------------------------------------
!  Having completed the loop for T,rho (in kit < kmax iterations) or    
!  discovered that the solution will not converge (making kit >= kmax   
!  depending on how it fails) the subroutine returns c  for the next set
!  of T, rho and Ye.  If ider >1 the subroutine calculates several momen
!  of the distribution, including the average A and Z and the binding en
!  ---------------------------------------------------------------------
      If(ider.ge.1) Then                                                
          dkT2=1.d0/(8.6174d-2*t9*t9)                                   
          atst=0.d0                                                     
          ztst=0.d0                                                     
          ytst=0.d0                                                     
          benuc=0.d0                                                    
          dt9=1.d0/t9                                                   
          DO 40 m=1,nmax                                                
            iz=intz(m)                                                  
            zpm=zp(m)                                                   
            zam=za(m)                                                   
            bim=bi(m)                                                   
            ym=dexp(dlgalp(m)+dlgyp*zpm+dlgyn*zn(m)+h(iz))              
!           y(m)=ym                                                     
            xm=zam*ym                                                   
!           x(m)=xm                                                     
            zpym=zpm*ym                                                 
!           at(iz)=at(iz)+xm                                            
!           zt(iz)=zt(iz)+zpym                                          
!                                                                       
!--abar, zbar, benuc, are the average baryon number, proton number      
!  number, and binding energy per nucleon (in Mev/nuc).                 
!                                                                       
            ytst=ytst+ym                                                
            atst=atst+xm                                                
            ztst=ztst+zpym                                              
            benuc=benuc+bim*ym                                          
   40     CONTINUE                                                      
!         Do 41 i=1,32                                                  
!           fye=zt(i)/at(i)                                             
!           write(15,2050) i,fye,zt(i),at(i)                            
!  41     Continue                                                      
          abar=atst/ytst                                                
          zbar=ztst/ytst                                                
!                                                                       
!--convert from binding energy per nucleon to ergs per gram             
!--with ye correction                                                   
          ubind=-9.616221376d17*(benuc+ye*0.783d0)                      
!         write(10,2000)(y(i),i=1,nmax)                                 
!         write(10,2100)atst,ztst                                       
      Endif                                                             
      Return                                                            
  999 Format(a5,6(1x,1pd10.3))                                          
 1000 FORMAT(3x,d9.2,6x,d9.2)                                           
 1100 FORMAT(5x,d11.4,8x,d10.3,5x,i3,8x,i3)                             
 1200 FORMAT(4x,d13.6,5x,d13.6,4x,d11.4)                                
 1300 FORMAT(5x,d9.2,5x,i3,4x,f4.1,4x,f4.1,4x,f4.1,4x,f4.2,4x,f4.2,     &
      4x,f4.2)                                                          
 2000 FORMAT (4(1pd15.7))                                               
 2050 Format(i2,3(1x,1pd10.3))                                          
 2100 FORMAT(1x,'Mass Conserv=',1pd15.7,' and Ye=',1pd15.7)             
 2200 Format(1x,6d11.3)                                                 
 2300 Format(a5,1x,d9.2,1x,a5,d9.2,1x,a5,d9.2,1x,a5,1x,d9.2)            
 3000 FORMAT(1x,'Convergence failure; Det=0 in step',I3,' at T9=',      &
      1pd9.2,' and density=',1pd9.2)                                    
 3100 FORMAT(1x,'Conv successful for T9=',1pd11.4,                      &
      ' rho=',1pd11.4,' ye=',1pd11.4)                                   
 3200 FORMAT(1x,'Convergence failure; does not converge in',I3,         &
      ' steps at T9=',1pd9.2,' and density=',1pd9.2)                    
 4000 FORMAT(1x,'k=',i3,' yp=',1pd16.8,' yn=',1pd16.8)                  
 4100 Format(1x,' delyp=',1pd13.6,' delyn=',1pd13.6,' det=',1pd13.6)    
 4200 Format(1x,' dfp=',1pd11.4,' dfn=',1pd11.4,' dgp=',1pd11.4,' dgn=',&
      1pd11.4)                                                         
 4250 Format(1x,' dfp2=',1pd11.4,' dfn2=',1pd11.4,' dgp2=',1pd11.4,     &
      ' dgn2=',1pd11.4)                                                
 4300 format(1x,'f=',1pd11.4,'g=',1pd11.4)                              
 4400 Format(1x,f4.1,1x,f4.1,1x,f7.4,1x,1pd10.3,1x,'scr')               
      END                                           
!                                                                       
      subroutine nucdata 
      implicit double precision (a-h,o-z) 
      double precision me 
      parameter (nyab= 16,nyabmx= 16) 
      dimension za(nyabmx),zp(nyabmx),zn(nyabmx),bi(nyabmx),angm(nyabmx) 
      dimension gt(24),dlg(24,nyabmx),it9i(24),t9i(24) 
      character*5 inam(nyabmx),name 
      common /zanp/ za,zp,zn,pnp,bi,angm 
      common /name/ inam 
      common /partf/ dlg,t9i 
      OPEN(8,FILE='netwinv4',STATUS='old') 
      read( 8,1000) name 
      read( 8,1010) (it9i(i),i=1,24) 
      DO 100 i=1,24 
        t9i(i)=it9i(i)*0.01 
  100 END DO 
      t9i(24)=t9i(24)*10. 
      DO 110 n=1,nyab 
        read( 8,1000) inam(n) 
  110 END DO 
!-----------------------------------------------------------------------
!  This subroutine reads, from the file netwinv3, the nuclear data which
!  be needed for later calcultions.  This data includes the atomic numbe
!  the number of protons and neutrons, and the binding energy (calculate
!  from the tabulated mass excess.  Also the tabulation of the  partitio
!  function, g, are read in for later interpolation. Once the set of nuc
!  data is read in, it is assigned to the proper nuclei.                
!-----------------------------------------------------------------------
      DO 120 l=1,nyab 
        read( 8,2000) name,a,na,nb,sp,me 
        read( 8,2001) (gt(m),m=1,24) 
        za(l)=a 
        zp(l)=dfloat(na) 
        zn(l)=dfloat(nb) 
        bi(l)=8.07144d0*zn(l)+7.28899d0*zp(l)-me 
        angm(l)=2.d0*sp+1.d0 
        inam(l)=name 
        DO 130 m=1,24 
          dlg(m,l)=dlog(gt(m)) 
  130   CONTINUE 
  120 END DO 
      bi(1)=0.d0 
      bi(2)=0.d0 
!     write(15,3000) (inam(i),za(i),zp(i),zn(i),bi(i),i=1,nyab)         
  999 Format(6(1x,1pd12.5)) 
 1000 format(a5) 
 1010 format(24i3) 
 2000 format(a5,f12.3,2i4,f6.1,f10.3) 
 2001 format(8f9.2) 
 3000 format(1x,a5,3f15.2,f15.3) 
      RETURN 
      END                                           
!                                                                       
      subroutine alphcalc (nmax,t9,rho) 
      implicit double precision (a-h,o-z) 
      parameter (nyabmx= 16,pi=3.14159,hbr=6.58217D-22,amu=1.036435E-18) 
      parameter (nyab= 16,bok=8.6174D-02,avn=6.02205D+23,c1=-77.7675d0) 
      dimension za(nyabmx),zp(nyabmx),zn(nyabmx),bi(nyabmx),angm(nyabmx) 
      dimension dlg(24,nyabmx),t9i(24) 
      dimension dlgalp(nyabmx) 
      character*5 inam(nyabmx) 
      common /zanp/ za,zp,zn,pnp,bi,angm 
      common /alph/ dlgalp 
      common /name/ inam 
      common /partf/ dlg,t9i 
!-----------------------------------------------------------------------
!  The first step in calculating alph is to interpolate the correct g,  
!  the partition function, for the temperature.                         
!-----------------------------------------------------------------------
      do j=1,24 
         i=j 
         if(t9i(j).gt.t9)go to 225 
      enddo 
  225 continue 
!                                                                       
      bkt=bok*t9 
      dbkt=1.d0/bkt 
!     tmp=1.5d0*dlog(((avn*rho)**(2.d0/3.d0)*2.d0*pi*hbr**2)/(bkt*amu)) 
      tmp=dlog(rho*avn)+c1-1.5d0*dlog(T9) 
!     If(t9.eq.5.) Then                                                 
!         write(15,*)rho,avn,hbr,pi,bkt,tmp                             
!     Endif                                                             
      dlgalp(1)=0.d0 
      dlgalp(2)=0.d0 
      ii=i 
      IF(ii.eq.1) THEN 
          DO 230 i=3,nmax 
            zat=za(i) 
            tmpa=dlog(angm(i)*zat*dsqrt(zat))+dlg(1,i)-.6931472d0*zat 
            dlgalp(i)=tmpa+(zat-1.d0)*tmp+bi(i)*dbkt 
  230     CONTINUE 
      ELSEIF(t9.ge.t9i(24)) THEN 
          DO 240 i=3,nmax 
            zat=za(i) 
            tmpa=dlog(angm(i)*zat*dsqrt(zat))+dlg(24,i)-.6931472d0*zat 
            dlgalp(i)=tmpa+(zat-1.d0)*tmp+bi(i)*dbkt 
  240     CONTINUE 
      ELSE 
          dt9i=t9i(ii)-t9i(ii-1) 
          ddt9i=1.d0/dt9i 
          dt9=t9-t9i(ii-1) 
          DO 250 i=3,nmax 
            grad=(dlg(ii,i)-dlg(ii-1,i))*ddt9i 
            gg=dexp(dlg(ii-1,i)+dt9*grad) 
            zat=za(i) 
            tmpa=dlog(angm(i)*gg*zat*dsqrt(zat))-.6931472d0*zat 
            dlgalp(i)=tmpa+(zat-1.d0)*tmp+bi(i)*dbkt 
  250     CONTINUE 
      ENDIF 
      RETURN 
  500 FORMAT (10f7.3) 
 1000 format(2a5) 
 1010 format(24i3) 
      END                                           
!                                                                       
      SUBROUTINE nsescreen(T9,rho,ye,h,dhdp,dhdn,dhdt9,iscrn) 
      implicit double precision (a-h,o-z) 
      dimension h(37),dhdp(37),dhdn(37),dhdt9(37) 
!                                                                       
      third=1.d0/3.d0 
!                                                                       
      GNP=2.27472d-4*(rho*ye)**third/T9 
      dGNPdT9=-GNP/t9 
      Htot=0.d0 
      dHtdp=0.d0 
      dHtdn=0.d0 
      dHtdT=0.d0 
      h(1)=0.d0 
      dhdp(1)=0.d0 
      dhdn(1)=0.d0 
      dhdT9(1)=0.d0 
      Do 100 j=1,36 
        z1=1.d0 
        z2=float(j) 
        a1=1.d0 
        a2=2.*z2 
        If(z2.eq.1.d0) a2=1.d0 
        ZZ=Z1*Z2 
        Z13=Z1**third 
        Z23=Z2**third 
        ZZm=2.d0*ZZ/(z13+z23) 
!-----STRONG SCREENING BY ITOH ET AL.(1990)                             
        GNZ=ZZm*GNP 
        dGNZdT9=ZZm*dGNPdT9 
!       GNZ=GNP*ZZ*2.d0/(Z13+Z23)                                       
        EM2=A1*A2*2.d0/(A1+A2) 
        TAU=3.3722d0*(EM2*ZZ*ZZ/T9)**third 
        dTAUdT9=TAU/(-3.d0*T9) 
        GT=GNZ/TAU 
        GT3=3.d0*GT 
        GT33=GT3*GT3*GT3 
        GT36=GT33*GT33 
        GT312=GT36*GT36 
        FNUM=.0455d0*GT3+.348d0*GT33+9.49d0*GT36-.123d0*GT312+          &
     &    .101d0*GT312*GT3                                              
        FDEN=1.d0+100.d0*GT33*GT3+.267*GT312 
        F90=FNUM/FDEN 
        SCR=1.25d0*GNZ-TAU*F90 
        SCRO=GNZ*(1.25d0-0.855d0*GT) 
        H0=SCR 
        IF (iscrn.eq.2) Then 
            H0=1.25d0*GNZ 
            dH0dT9=1.25d0*dGNZdT9 
            DHDYP=0.d0 
            DHDYN=0.d0 
!           write(*,*) 'H',h0,dh0dt9                                    
        Endif 
        FSR=GNZ*(1.25d0-1.71d0*GT)*0.333333333333333333d0 
        FST=-GNZ*(1.25d0-1.425d0*GT) 
   30   Continue 
!-----                                                                  
!  Add succeeding screening factors                                     
        Htot=Htot+h0 
        dHtdp=dHtdp+dhdyp 
        dHtdn=dHtdn+dhdyn 
        dHtdT9=dHtdT9+dH0dT9 
        h(j+1)=Htot 
        dhdp(j+1)=dHtdp 
        dhdn(j+1)=dHtdp 
        dhdT9(j+1)=dHtdT9 
        scrfct=exp(h(j+1)) 
!       If(T9.eq.Tfnl.and.rho.eq.rhofnl) Then                           
!           write(15,999) 'scrn',z1,z2,h0,dh0dp,dh0dn                   
!           tstscr=exp(h0)                                              
!           write(15,*) 'TAU',tau,gnz                                   
!           write(15,999) 'total',h(j+1),dhdp(j+1),dhdn(j+1)            
!           write(14,999) 'scrn',z1,z2,h0,tstscr,scrfct                 
!       Endif                                                           
  100 Continue 
      RETURN 
  999 Format(a5,6(1x,1pd10.3)) 
      END                                           
!                                                                       
      subroutine newnorm(x,za,n1,n,zp,zn,ye) 
      implicit double precision (a-h,o-z) 
      dimension x(n),za(n),zp(n),zn(n) 
      a=0.0d0 
      b=0.0d0 
      c=0.0d0 
      d=0.0d0 
!-----------------------------------------------------c                 
!   summation over loop n1 covers read-in abundances  c                 
!   if mass fractions ommit za(i)                     c                 
      Do 10 i=1,n1 
        a=a+zn(i)*x(i) 
        b=b+zp(i)*x(i) 
        c=c+zn(i)*zp(i)*x(i)/za(i) 
        d=d+(zp(i)**2)*x(i)/za(i) 
!       write(6,*) x(i),za(i),zp(i),zn(i),a,b                           
   10 Continue 
!     write(15,*) 'newnorm',a,b,c,d                                     
!-----------------------------------------------------c                 
      n2=n1+1 
      s2=0.d0 
      if(n2.gt.n) goto 18 
!-----------------------------------------------------                  
!   takes solar abundances from n2 to n                                 
!-----------------------------------------------------                  
      do 15 i=n2,n 
        a=a+zn(i)*x(i) 
        b=b+zp(i)*x(i) 
        c=c+zn(i)*zp(i)*x(i)/za(i) 
        d=d+zp(i)**2*x(i)/za(i) 
   15 Continue 
   18 beta=(ye*a-c)/(a*d-b*c) 
      alph=(1-beta*b)/a 
!-----------------------------------------------------c                 
      Do 20 i=1,n1 
        x(i)=x(i)*(alph*zn(i)+beta*zp(i))/za(i) 
   20 Continue 
!-----------------------------------------------------c                 
      return 
      END                                           
                                                                        
!********************************************************               
!                                                                       
!     Nadyozhin eos, obtained through Stan Woosley                      
!                                                                       
!*******************************************************                
!                                                                       
      subroutine nados(tt,dd,zbar,abar,pel,eel,sel,                     &
     &                 ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)         
      implicit real*8 (a-h,o-z) 
!..                                                                     
!..                                                                     
      double precision tt,dd,zbar,abar,pel,eel,sel 
!..                                                                     
!..communicate                                                          
      common/arg/t,den,psi 
      common/iarg/lst,kentr,kpar,jurs,jkk 
      common/nz/nz 
      common/az/as,zs,scn 
      common/result/p,e,s,sk,pt,et,st 
      common/resel/pe,ee,se,sek,hpr 
      common/str/ppl,epl,spl,cp,gam,da,dpe,dse,dsp,beta 
      common/nzr/nzr 
!..                                                                     
!..t in 10**9 den in 10**7                                              
!     t   = tt * 1.0e-9                                                 
      t   = tt 
!     den = dd * 1.0d-7                                                 
      den = dd 
      as  = abar 
      zs  = zbar 
!..      scn = 10.063379                                                
      scn = 2.5 * log(abar) 
!                                                                       
!..get temp and density derivatives; get entropy; get number of pairs   
!-- we need the works                                                   
      nz    = 0 
      jurs  = 0 
      lst   = 2 
      kentr = 1 
      kpar  = 1 
!..                                                                     
      call epeos 
!..                                                                     
!..return arguments                                                     
!     pel  = pe * 1.0d24                                                
!     eel  = ee * 1.0e17                                                
!     sel  = se * 1.0e8                                                 
!     ptot = p  * 1.0d24                                                
!     etot = e * 1.0e17                                                 
!     stot = s * 1.0e8                                                  
      pel  = pe 
      eel  = ee 
      sel  = se 
      ptot = p 
      etot = e 
      stot = s 
      dpt  = pt 
      det  = et 
      dpd  = ppl 
      ded  = epl 
      gamm = gam 
!--somehow, eta is returned negative in the perfect gas asymptotic      
      eta  = dmax1(psi,0.d0) 
!..                                                                     
!..      write(6,2) t,den,lst,kentr,nz,nzr                              
!..      write(6,4) p,e,s,sk                                            
!..      write(6,5) pt,et,st,ppl                                        
!..      write(6,6) epl,spl,pe,ee                                       
!..      write(6,7) se,sek,hpr,gam                                      
!..      write(6,8) da,dpe,dse,dsp                                      
!..      write(6,9) psi,beta                                            
!..                                                                     
!..2     format(3x,'t=',d10.3,4x,'den=',d10.3,4x,'lst=',i1,4x,          
!..     * 'kentr=',i1/3x,'nz=',i3,'  nzr=',i3/1x)                       
!..4     format(1x,'p  =',d12.5,'  e   =',d12.5,'  s  =',d12.5,         
!..     * '  sk =',d12.5)                                               
!..5     format(1x,'pt =',d12.5,'  et  =',d12.5,'  st =',d12.5,         
!..     * '  ppl=',d12.5)                                               
!..6     format(1x,'epl=',d12.5,'  spl =',d12.5,'  pe =',d12.5,         
!..     * '  ee =',d12.5)                                               
!..7     format(1x,'se =',d12.5,'  sek =',d12.5,'  hpr=',d12.5,         
!..     * '  gam=',d12.5)                                               
!..8     format(1x,'da =',d12.5,'  dpe =',d12.5,'  dse=',d12.5,         
!..     * '  dsp=',d12.5)                                               
!..9     format(1x,'psi=',d12.5,'  beta=',d12.5)                        
!..                                                                     
!..                                                                     
! *** an example of calculation of half-integer fermi-dirac functions   
! *************** results are in common/fdf/ -- f-d functions and deriva
!..      psi=3.d0                                                       
!..      call fd12f                                                     
      return 
      END                                           
!..                                                                     
!..                                                                     
!..                                                                     
!..                                                                     
!..                                                                     
      subroutine epeos 
!          ***  version 1.1 santa cruz, august 2, 1992  ***             
!***********************************************************************
!  *** equation of state for completely ionized matter                  
!  *** electron & positron component --- fermi-dirac statistics using   
!                       various asymptotic expansions where possible    
!  *** ion component --- a perfect gas approximation                    
!  *** black-body radiation                                             
!***********************************************************************
!                            references                                 
!   1. nadyozhin d.k. 1974, "naucnye informatsii", nos. 32, 33          
!                            (ucsc science library: qb 1 a4)            
!***********************************************************************
       implicit real*8 (a-h,o-z) 
!                                                                       
!  *** the arguments                                                    
      common/arg/t,den,psi 
      equivalence(den,pl) 
      common/iarg/lst,kentr,kpar,jurs,jkk 
!***********************************************************************
!                                                                       
!  ***   t --- temperature in 10**9 k                                   
!  *** den --- density in 10**7 g/ccm                                   
!  *** psi --- parameter of degeneracy. works as an argument only when  
!              one enters entry fd12f to get fermi-dirac functions,     
!              otherwise it is calculated as a function of t and den.   
!  *** lst, kentr, kpar --- the keys aimed to make the calculations     
!               faster where possible (default: lst, kentr, kpar = 0)   
!  *** lst=0 --- epeos calculates only thermodynamic functions p,e,s    
!      lst=1 --- the first temperature-derivatives are calculated       
!                in addition to the thermodynamic functions             
!      lst=2 --- extends calculations to get density-derivatives as well
!  *** kentr=0 --- turns the calculation of entropy off and suspends    
!            the calculation of psi in a perfect gas asymptotics (nz=1).
!      kentr=1 --- turns the calculation of entropy on.                 
!  *** kpar=0 --- when in relativistic asymptotics (nz=4), turns off the
!        calculation of total number of pair (hpr),(kpar=1 --- turns on)
!      kpar is inactive for other asymptotics.                          
!                                                                       
!  *** jkk --- the current number of mesh point, is inactive in this    
!        version of epeos. however, it appears in epeos error messages  
!        and, if defined, may be useful for locating of errors in       
!        a program calling epeos.                                       
!***********************************************************************
      common/nz/nz 
!***********************************************************************
!  *** nz --- specifies the operational mode (default: nz=0)            
!      nz=0 --- calculations with the overall search for                
!               the appropriate working region                          
!      for 0<nz<6 epeos works within one of five following modes        
!      independent of the values of temperature and density specified   
!      nz=1 --- perfect gas approximation with the first order          
!               corrections for degeneracy and pairs                    
!      nz=2 --- expansion over half-integer fermi--dirac functions      
!      nz=3 --- chandrasekhar's expansion for degenerate gas            
!      nz=4 --- relativistic asymptotics                                
!      nz=5 --- quadratures taken with the gauss method                 
!***********************************************************************
      common/az/as,zs,scn 
!***********************************************************************
!  ***  as --- mean mass number, zs --- mean nuclear charge             
!          emue=as/zs --- mean electron molecular weight: the total     
!          number of "atomic" electrons in a unit volume, nea,          
!          amounts to (density)/(mu*emue), mu is atomic mass unit.      
!          for a mixture:   as=1/sum{y_i},  zs=as*sum{z_i*y_i}, where   
!          y_i=x_i/a_i,  and  x_i  being a mass fraction of 'i' species.
!  *** scn --- additive entropy constant for the ion component          
!          for a gas of nuclei with mass number a, and spin i           
!          scn=ln[(2i+1)*a**2.5], for iron-56, scn=2.5*ln(56)=10.063379.
!          for a mixture:   scn=as*sum{y_i*ln[(2i_i+1)*(a_i)**2.5}.     
!  *** jurs --- if jurs=0 then common-block  /az/ is to be preliminary  
!            filled with all necessary information. (default values     
!            of as,zs,scn are specified in block data for pure iron-56).
! if jurs=1 (default), epeos applies to subroutine chemic for as,zs,scn 
!***********************************************************************
!                                                                       
!                         diagnostics                                   
!***********************************************************************
! epeos opens file 'epeos.mes', writes the epeos version in it, analyses
! the values of arguments in common-blocks /arg/, /iarg/, /nz/, /az/ and
! then writes additional information in 'epeos.mes' if something wrong  
! with the arguments: in particular, epeos stops when  t = or < 0.      
! in case  den < 0 , epeos sends a warning in 'epeos.mes' and continues 
! to calculate with den=0 (black body radiation only).                  
!***********************************************************************
!                                                                       
!  *** the results of calculations                                      
      common/result/p,e,s,sk,pt,et,st 
      common/str/ppl,epl,spl,cp,gam,da,dpe,dse,dsp,beta 
      common/resel/pe,ee,se,sek,hpr 
      common/nzr/nzr 
!***********************************************************************
!  *** p --- total pressure in units 10**24 erg/ccm                     
!  *** e --- total specific energy in units 10**17 erg/g                
!  *** s --- total entropy in units  10**8 erg/(g k)                    
!  *** sk --- total dimensionless entropy per nucleon: sk=s*mu/kb       
!             mu -- atomic mass unit, kb -- boltzmann's constant        
!  *** pt,et,st --- temperature derivatives of p,e,s at constant density
!  *** ppl,spl --- density derivatives of p and s at constant temperatur
!  *** epl --- density derivatives of e multiplied by density den       
!  *** pe,ee,se,sek ----  pressure, specific energy and entropy         
!                         of the electron-positron gas component        
!  *** hpr --- the total number of the electron-positron pairs          
!              per "atomic" electron (hpr=mu*emue*np/den), where        
!              'np' is the number of pairs per unit volume.             
!  *** cp --- specific heat at constant pressure                        
!  *** gam --- adiabatic index = {d log(p)/d log(den)} at s=const       
!  *** da --- logarithmic adiabatic temperature gradient                
!  *** dpe = (den/p)(epl+t*pt/den)-1 -- thermodynamic identity: dpe=0   
!  *** dse = t*st/et-1 -- thermodynamic identity: dse=0                 
!  *** dsp = -spl*(den/pt)-1 -- thermodynamic identity: dsp=0           
!  *** beta --- ratio of black-body radiation pressure to total pressure
!  *** nzr --- identificator of working region on t-den plane when nz=0 
!***********************************************************************
      common/fdf/f12,f32,f52,f72,f12s,f32s,f52s,f72s 
!***********************************************************************
!  *** f12,f32,f52,f72 --- half-integer fermi-dirac functions           
!  *** f12s,f32s,f52s,f72s --- the first derivatives of f-d functions   
!  *** psi (in common/arg/t,den,psi) is the argument of f-d functions   
!      there exists special entry to get these functions separately --  
!      use command call fd12f after specifying psi in common/arg/       
!***********************************************************************
      dimension cu(62),ck(3),a(17),c(8),b(22),b1(4),b2(4),b3(4),b4(4),  &
     &b5(4),c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),d1(4),d2(4),d3(4),d4(4),&
     &d5(4),d6(4),d(4),a1(4),a2(4),a3(4),a4(4),df1(4),df2(4)            
      dimension uio(5),ui1(5),ui2(5),cio(5),ci1(5),ci2(5),aio(5),ai1(5) 
      dimension ai2(5),xxi(5),aai(5),cci(5),bbi(5),wk1(5),wk2(5),wk3(5) 
      dimension uac(95),wk4(5),wk5(5),wk6(5) 
      dimension cpp(5),abc(85),ado(5),ad1(5),ad2(5),bdo(5),bd1(5),fgs(8) 
      dimension bd2(5),cdo(5),cd1(5),cd2(5),gdo(5),gd1(5),gd2(5) 
      dimension ggsi(5),zzi(5),vvi(5),hhi(5),ggi(5) 
      dimension asp(3),bsp(3),csp(3),gsp(3),aspa(3),bspa(3),cspa(3),gspa&
     &(3),abcg(24),wk7(5)                                               
      equivalence(psi,pc) 
      equivalence (uio(1),uac(1)),(ui1(1),uac(6)),(ui2(1),uac(11)),     &
     &(cio(1),uac(16)),(ci1(1),uac(21)),(ci2(1),uac(26)),(aio(1),uac(31)&
     &),(ai1(1),uac(36)),(ai2(1),uac(41)),(xxi(1),uac(46)),(aai(1),uac(5&
     &1)),(cci(1),uac(56)),(bbi(1),uac(61)),(wk1(1),uac(66)),(wk2(1),uac&
     &(71)),(wk3(1),uac(76)),(wk4(1),uac(81)),(wk5(1),uac(86)),(wk6(1),u&
     &ac(91))                                                           
      equivalence(abc(1),ado(1)),(abc(6),ad1(1)),(abc(11),ad2(1)),      &
     &(abc(16),bdo(1)),(abc(21),bd1(1)),(abc(26),bd2(1)),(abc(31),cdo(5)&
     &),(abc(36),cd1(1)),(abc(41),cd2(1)),(abc(46),gdo(1)),(abc(51),gd1(&
     &1)),(abc(56),gd2(1)),(abc(61),ggsi(1)),(abc(66),zzi(1)),          &
     &(abc(71),vvi(1)),(abc(76),hhi(1)),(abc(81),ggi(1)),(abcg(1),asp(1)&
     &),(abcg(4),bsp(1)),(abcg(7),csp(1)),(abcg(10),gsp(1)),(abcg(13),as&
     &pa(1)),(abcg(16),bspa(1)),(abcg(19),cspa(1)),(abcg(22),gspa(1))   
      data eit/1.d-6/,dst/1.d-3/,grcf/1.3d0/,grif/0.7d0/,tpar/0.3d0/,   &
     &tt1/0.8d0/,tt2/0.96d0/,tpar1/0.32d0/,ro1/0.12d0/,ro2/3.4884680d-2/&
     &ep2/0.7d0/,ps1/7.2427389d-5/,ro1l/-2.1202635d0/,ro2l/-3.3557075d0/&
     &xep2/0.3d0/,rotp/3.1220424d-10/,pss1/0.1662467d0/,                &
     & pss2/0.1881164d0/                                                
      data cu/5.93013d0,0.831434d0,1.5d0,1.d0,2.8125d0,5.5957031d0,     &
     &3.75d0,5.15625d0,1.25d0,0.9375d0,1.8652344d0,5.25d0,1.265625d1,   &
     &0.5625d0,2.d0,0.5d0,3.d0,2.25d0,0.75d0,2.984375d0,1.3125d1,       &
     &1.6875d1,3.5d0,2.5d0,1.6787109d1,1.40625d0,2.5d0,0.25d0,          &
     &5.6568542d0,5.75d0,1.453125d1,4.63125d1,5.8125d1,1.35d1,          &
     &1.013429d-3,1.33333333333d0,4.d0,2.8613184d-2,0.94051952d-1,      &
     &0.71428571d0,0.39650038026d-1,0.21875d0,0.66666666667d0,0.6d0,    &
     &1.05d0,0.18007375d0,0.292178d0,0.3601475d0,0.33333333333d0,       &
     &5.d0,8.d0,2.66666666667d0,24.d0,15.d0,6.d0,9.d0,1.2d1,            &
     &1.02677135171d+1,4.9305117064d0,4.80195683116d-1,                 &
     &8.09755744169d-2,7.99196885645d-1/                                
      data a1/6.78196814d-1,6.78183752d-1,6.94779549d-1,7.60042563d-1/ 
      data a2/5.37738527d-1,5.37346659d-1,4.87559266d-1,4.09243650d-1/ 
      data a3/1.73981666d-1,1.70062988d-1,2.19850381d-1,0.251176627d0/ 
      data a4/2.38231503d-2,1.07608902d-2,-5.83490747d-3,-.0100117403d0/ 
      data df1/3.47963332d-1,3.40125976d-1,4.39700762d-1,0.502353254d0/ 
      data df2/7.14694509d-2,3.22826706d-2,-1.75047241d-2,              &
     &         -3.00352209d-2/                                          
      data  b1/1.15292705d0,1.15292656d0,1.14670313d0,1.08551906d0/ 
      data  b2/1.01729522d0,1.01727563d0,1.04216932d0,1.14006384d0/ 
      data  b3/4.03303895d-1,4.03009994d-1,3.65669449d-1,3.06932737d-1/ 
      data b4/8.69908330d-2,8.50314940d-2,1.09925190d-1,1.25588313d-1/ 
      data b5/8.93368136d-3,4.03533382d-3,-2.18809030d-3,-3.75440261d-3/ 
      data c1/3.08286343d0,3.08286341d0,3.08597512d0,3.16245521d0/ 
      data c2/2.88231761d0,2.88231639d0,2.86675783d0,2.71379764d0/ 
      data c3/1.27161903d0,1.27159453d0,1.30271165d0,1.42507981d0/ 
      data c4/3.36086579d-1,3.35841662d-1,3.04724541d-1,2.55777281d-1/ 
      data c5/5.43692706d-2,5.31446837d-2,6.87032441d-2,7.84926959d-2/ 
      data c6/4.46684068d-3,2.01766691d-3,-1.09404515d-3,-1.87720131d-3/ 
      data d1/1.11841602d1,1.11841602d1,1.11823451d1,1.10708116d1/ 
      data d2/1.07900220d1,1.07900219d1,1.08009129d1,1.10685932d1/ 
      data d3/5.04405582d0,5.04405368d0,5.01682620d0,4.74914588d0/ 
      data d4/1.48355553d0,1.48352696d0,1.51983026d0,1.66259311d0/ 
      data d5/2.94075757d-1,2.93861454d-1,2.66633974d-1,2.23805121d-1/ 
      data d6/3.80584894d-2,3.72012786d-2,4.80922708d-2,5.49448872d-2/ 
      data d/2.60565706d-3,1.17697237d-3,-6.38193005d-4,-1.09503410d-3/ 
      data a/-3.53553400d-1,1.92450000d-1,-1.25000000d-1,8.94427000d-2  &
     &      ,-1.76777000d-1,6.41500000d-2,-3.12500000d-2,1.78885000d-2  &
     &      ,-8.83883000d-2,2.13833000d-2,-7.81250000d-3,3.57771000d-3  &
     &      ,1.16317000d1,-4.41942000d-2,7.12778000d-3,-1.95313000d-3   &
     &      ,7.15541750d-4/                                             
      data b/6.666667d-1,8.22467d-1,7.102746d-1,6.467679d0              &
     &      ,4.00000000d-1,2.46740000d0,-7.10275000d-1,-2.77186000d0    &
     &      ,2.85714000d-1,4.11234000d0,3.55137000d0,2.77186000d0       &
     &      ,2.22222200d-1,5.75727000d0,2.48596000d1,-6.4676800d0       &
     &      ,7.79429075d-3,-4.94746507d-2,1.94857269d-2,3.56600189d-2   &
     &      ,-1.73161277d-1,3.41000220d-2/                              
      data c/-7.07106800d-1,5.77350000d-1,-5.00000000d-1,4.47213500d-1  &
     &      ,1.00000015d0,-4.11233500d-1,-1.77568650d0,-2.91045555d1/   
      data ck/8.86227d-1,1.32934d0,3.32335d0/ 
      data uio/0.43139881d0,1.7597537d0,4.1044654d0,7.7467038d0,        &
     &         1.3457678d1/                                             
      data ui1/0.81763176d0,2.4723339d0,5.1160061d0,9.0441465d0,        &
     &         1.5049882d1/                                             
      data ui2/1.2558461d0,3.2070406d0,6.1239082d0,1.0316126d1,         &
     &         1.6597079d1/                                             
      data cio/0.37045057d0,0.41258437d0,9.777982d-2,5.3734153d-3,      &
     &         3.8746281d-5/                                            
      data ci1/0.39603109d0,0.69468795d0,0.2232276d0,1.5262934d-2,      &
     &         1.3081939d-4/                                            
      data ci2/0.76934619d0,1.7891437d0,0.70754974d0,5.6755672d-2,      &
     &         5.557148d-4/                                             
      data aio/0.64959979d0,0.17208724d0,0.016498837d0,4.321647d-4,     &
     &         1.4302261d-6/                                            
      data ai1/0.44147594d0,8.4387677d-2,5.9999383d-3,1.180802d-4,      &
     &         2.9101763d-7/                                            
      data ai2/0.28483475d0,4.0476222d-2,2.1898807d-3,3.3095078d-5,     &
     &         6.194128d-8/                                             
      data xxi/7.265351d-2,0.2694608d0,0.5331220d0,0.7868801d0,         &
     &         0.9569313d0/                                             
      data aai/3.818735d-2,0.1256732d0,0.1986308d0,0.1976334d0,         &
     &         0.1065420d0/                                             
      data cci/0.26356032d0,1.4134031d0,3.5964258d0,7.0858100d0,        &
     &         1.2640801d1/                                             
      data bbi/0.29505869d0,0.32064856d0,7.391557d-2,3.6087389d-3,      &
     &         2.3369894d-5/                                            
      data pc1/0.5d0/,pc2/0.7d0/ 
      data cpp/5.d0,1.d1,1.5d1,2.5d1,2.5d0/ 
      data fgs/0.571428571d0,0.33333333d0,0.2272727d0,0.168269d0,       &
     & 0.142857143d0,5.5555555d-2,2.840909d-2,1.68269d-2/               
      data (abc(k),k=1,76)/7.72885519d1,                                &
     &  1.42792768d2, 4.30552910d1, 2.43440537d0,                       &
     &  1.75674547d-2, 9.99400362d1, 2.73430265d2, 1.00130386d2,        &
     &  6.91871969d0, 5.93132645d-2, 2.30460043d2, 7.56122303d2,        &
     &  3.19543255d2, 2.57313963d1, 2.51960145d-1,-2.35425685d1,        &
     & -4.38697759d1,-1.32985534d1,-7.52438243d-1,-5.42994019d-3,       &
     & -3.05287674d1,-8.42357074d1,-3.09412790d1,-2.13850223d0,         &
     & -1.83331906d-2,-7.06062732d1,-2.33317365d2,-9.87584116d1,        &
     & -7.95332909d0,-7.78785901d-2, 1.42401164d-1, 4.12778210d-1,      &
     &  1.52786427d-1, 8.84665279d-3, 6.38815164d-5, 2.18702630d-1,     &
     &  8.82651141d-1, 3.60865258d-1, 2.51545288d-2, 2.15684504d-4,     &
     &  5.87304073d-1, 2.59226969d0, 1.15817403d0, 9.35640728d-2,       &
     &  9.16218624d-4, 2.94914091d-1, 5.29893251d-1, 1.56942521d-1,     &
     &  8.85295620d-3, 6.38816670d-5, 3.77889882d-1, 1.00545595d0,      &
     &  3.64435019d-1, 2.51594259d-2, 2.15684607d-4, 8.63109766d-1,     &
     &  2.76526224d0, 1.16235562d0, 9.35691781d-2, 9.16218718d-4,       &
     &  7.22774549d-1, 6.91700407d-1, 6.36940508d-1, 5.69038300d-1,     &
     &  5.14846835d-1, 9.63560320d-1, 2.11340310d0, 4.29642580d0,       &
     &  7.78581000d0, 1.33408010d1, 5.08574570d-2, 1.88622560d-1,       &
     &  3.73185400d-1, 5.50816070d-1, 6.69851910d-1, 2.89632880d-1/     
        data (abc(k),k=77,85)/                                          &
     &  4.66144392d-1, 1.53210873d-1, 1.00694874d-2, 8.53586810d-5,     &
     &  1.46896384d-2, 4.60107809d-2, 6.75861819d-2, 6.21820743d-2,     &
     &  3.16690579d-2/                                                  
          data nitm/30/ 
          data pi2/9.8696044011d0/ 
          data t5/1.3d1/,t4/1.5d1/,ro3/3.d1/,ro4/3.9d1/ 
          data nfil/1/ 
! ***                                                                   
        if(nfil.ne.1) go to 4000 
!..        open(unit=101,file='epeos.mes')                              
!..        write(6,5001)                                                
        nfil=0 
 4000 if((nz.ge.0).and.(nz.le.5)) go to 4004 
      write(6,5002) nz 
!..      print 4005,nz                                                  
!..      print 4006                                                     
      stop 
 4005 format('  illegal value of nz:  nz=',i5) 
 4006 format(' allowed values are nz=0,1,2,3,4,5  *stop in epeos*') 
 5001   format(10x,'epeos  ***  version 1.1 august 2, 1992  ***  epeos') 
 5002   format(20x,'epeos  ***  error in   nz   ***  epeos'/            &
     & 1x,'illegal value of nz:  nz=',i5,                               &
     & 1x,'! allowed values are nz=0,1,2,3,4,5 *stop*')                 
 5012   format(20x,'epeos  ***  error in   t    ***  epeos'/            &
     & 1x,1p,'temperature t must be positive: t=',d13.6,                &
     & 1x,'jkk=',i4,' *stop*')                                          
 5013   format('  temperature t must be positive: t=',                  &
     & 1p,d13.6,' jkk=',i4,' *stop in epeos*')                          
 5022   format(20x,'epeos  ***     warning!     ***  epeos'/            &
     & 1x,1p,'negative or zero density: den=',d13.6,' jkk=',i4/         &
     & 1x,'calculations are going on with den=0.')                      
 5023   format('  negative or zero density den=',                       &
     &  d13.6,' jkk=',i4,' *warning in epeos* ')                        
 5032   format(20x,'epeos  *** error in as or zs ***  epeos'/           &
     & 1x,1p,'as and zs must be positive: as=',d10.3,' zs=',d10.3,      &
     & 1x,' jkk=',i4,' *stop*')                                         
 5033   format(' as and zs must be positive: as=',1p,d10.3,             &
     & 1x,'zs=',d10.3,' jkk=',i4,'*stop in epeos*')                     
!                                                                       
 4004 continue 
      if(t.gt.0.d0) go to 4010 
      write(6,5012) t,jkk 
!..      print 5013,t,jkk                                               
      stop 
 4010 continue 
! *** preparation for calculations                                      
      if(pl.gt.0.d0) go to 102 
      write(6,5022) pl,jkk 
!..      print 5023,pl,jkk                                              
      fac=t*t*t 
      pt=3.025884d-2*fac/cu(17) 
!     pt=3.025884d-2*t**3/cu(17)                                        
      p=0.25d0*t*pt 
      ppl=0.d0 
      gam=cu(36) 
      da=0.25d0 
      go to 9 
  102 if(jurs.eq.1) stop 'tried a call to chemic' 
!..call chemic                                                          
! *** subroutine chemic calculates as,zs,scn                            
      if((as.gt.0.d0).and.(zs.gt.0.d0)) go to 4030 
      write(6,5032) as,zs,jkk 
!..      print 5033,as,zs,jkk                                           
      stop 
 4030 emue=as/zs 
      rg=cu(2)/emue 
      ki=1 
      hpr=0.d0 
   90 alf=cu(1)/t 
      al1=cu(4)/alf 
      plm=pl/emue 
      sqe=0.d0 
      ei=t*rg 
      pi=ei*pl 
      eg=cu(3)*ei 
  590  continue 
!                                                                       
! *** search for required working region                                
      if(nz.ne.0) go to(1,2,3,4,5),nz 
      if(ki.ne.1) go to  123 
      if(plm.le.ro3) go to 310 
      if(plm.ge.ro4) go to 4 
      if(t.le.t5) go to 550 
      if(t.ge.t4) go to 4 
! *** searching around the triangle                                     
          x=(ro4-ro3)/(t4-t5) 
          y=ro3+(t4-t)*x 
          if(y.gt.plm) go to 800 
          go to 4 
  310  if(t.le.t5) go to 123 
          if(t.ge.t4) go to 4 
! *** interpolation over t in region 45                                 
          t1=t5 
          t2=t4 
          nz2=4 
          nz1=5 
          go to 128 
  550  continue 
! *** interpolation over density for density < ro4                      
      nzp1=0 
      if(t.lt.tt2) nzp1=3 
  583  kin=ki 
          nz=nzp1 
          ki=6 
          go to 590 
  577  pn1=pe 
        en1=ee 
         sn1=se 
         psn1=psi 
          hprn=hpr 
         pnp1=ppl 
        enp1=epl 
          snp1=spl 
          pnt=pt 
         ent=et 
         snt=st 
         nzr1=nzr 
         nz=4 
         ki=7 
          go to 590 
  578  z1=ro4-ro3 
          x2=(plm-ro3)/z1 
          fac=x2*x2 
          x=fac*(3.d0-2.d0*x2) 
!         x=x2**2*(3.d0-2.d0*x2)                                        
          x1=1.d0-x 
          z1=z1*emue 
          x3=6.d0*x2*(1.d0-x2)/z1 
          p=pe 
         e=ee 
        s=se 
          pe=pe*x+pn1*x1 
          ee=ee*x+en1*x1 
          if(kentr.eq.0) go to 591 
          se=se*x+sn1*x1 
  591 if(lst.eq.0) go to 592 
          pt=pt*x+pnt*x1 
          et=et*x+ent*x1 
          ppl=ppl*x+pnp1*x1+(p-pn1)*x3 
          epl=epl*x+enp1*x1+(e-en1)*x3 
          if(kentr.eq.0) go to 592 
          st=st*x+snt*x1 
          spl=spl*x+snp1*x1+(s-sn1)*x3 
  592  psi=psi*x+psn1*x1 
          hpr=hpr*x+hprn*x1 
          nzr=10*nzr1+nzr 
          ki=kin 
          go to 134 
  800 continue 
! ***** the triangle                                                    
          nz2=4 
         nz1=0 
          t1=t5 
          t2=t4-(plm-ro3)/x 
          go to 128 
  123 if(ki.ne.4) go to 136 
      if(nz2.eq.5) go to 111 
      nzp1=5 
          go to 583 
  136 if(plm.lt.ps1) go to 110 
      if(t.lt.tt1) go to 111 
      kzz=2 
         go to 121 
  115 y=x*grcf 
      if(plm.gt.y) go to 3 
      kzz=3 
         z1=t 
      if(t.lt.tt2) go to 112 
      if(plm.lt.x) go to 5 
  113 kkz=0 
  116 go to 121 
  114 dl=(x-plm)/xz 
      x1=1.d0 
      if(dl.lt.0.d0) x1=-1.d0 
      x2=dl*x1 
      if(x2.gt.0.3d0) dl=0.3d0*x1 
      t=t*(1.d0-dl) 
      if(x2.gt.eit)go to 116 
      if(kkz.eq.1) go to 118 
      t2=t 
      t=z1 
  138 z2=plm 
         plm=plm/grcf 
         kkz=1 
      go to 116 
  118 t1=t 
         nz1=3 
         plm=z2 
         t=z1 
         nz2=5 
      go to 128 
  112 if(plm.gt.pss1) go to 117 
      t1=tt1 
         t2=tt2 
        nz1=0 
        nz2=5 
      go to 128 
  117 if(plm.gt.pss2) go to 113 
      t2=tt2 
      go to 138 
  110 if(t.lt.tpar) go to 111 
      if(t.gt.tt2) go to 5 
      sqe=exp(-alf) 
      if(t.gt.tt1) go to 119 
      if(plm.gt.ro1*sqe) go to 1 
      if(t.gt.tpar1) go to 119 
      if(plm.gt.rotp) go to 120 
      t1=tpar 
        t2=tpar1 
  122 nz1=1 
         nz2=5 
         sqe=0.d0 
      go to 128 
  119 if(plm.lt.ro2*sqe) go to 5 
  120 t1=log(plm) 
      t2=cu(1)/(ro2l-t1) 
         t1=cu(1)/(ro1l-t1) 
      go to 122 
  111 sq=t*sqrt(t) 
      y=2.095d-2*sq 
      x=y*grif 
      if(plm-x)44,44,51 
!                                                                       
! perfect gas with corrections for degeneracy and pairs (nz=1)          
    1 sq=t*sqrt(t) 
   44 nzr=1 
      qa=6.9712909d0*plm/sq 
      x=qa*al1 
      x2=qa-x*(cu(5)-cu(6)*al1) 
      pe=cu(4)+x2 
      epl=qa-x*(cu(10)+cu(11)*al1) 
      x8=cu(9)*al1 
      x3=x8*(cu(4)+(cu(14)*al1-cu(4))*al1) 
      ee=cu(4)+epl+x3 
      if(lst.eq.0) go to 6 
      pt=cu(4)-cu(16)*x2-x*(cu(5)-cu(15)*cu(6)*al1) 
      ppl=cu(4)+cu(15)*x2 
      et=cu(4)-cu(16)*epl+x8*(cu(15)+al1*(cu(18)*al1-cu(17))-qa*        &
     &   (cu(19)+cu(20)*al1))                                           
    6 if(t.lt.tpar) go to 8 
      if(sqe.eq.0.d0) sqe=exp(-alf) 
      fac=al1*sqe/plm 
      x4=0.268192d0*fac*fac 
!     x4=0.268192d0*(al1*sqe/plm)**2                                    
      x5=x4*al1 
      hpr=x5*(cu(4)+al1*(cu(7)+cu(8)*al1)) 
      pe=pe+hpr 
      x6=(x4+x5*(cu(12)+cu(13)*al1))/cu(3) 
      ee=ee+x6 
      if(lst.eq.0) go to 8 
      ppl=ppl-hpr 
      pt=pt+hpr*cu(15)*(cu(15)+alf)+x5*(cu(7)+cu(21)*al1)*al1 
      x8=x6*cu(15) 
      et=et+x8*(cu(3)+alf)+x5*(cu(23)+cu(22)*al1) 
      epl=epl-x8 
    8 if((kentr+ki).eq.1) go to 56 
      x7=cu(16)*(qa+x*(cu(5)-cu(25)*al1)) 
      psi=log(cu(29)*qa) 
      se=cu(24)-psi+x7+al1*(cu(7)+al1*(al1*cu(26)-cu(5))) 
      psi=psi+2.d0*x2+0.5d0*hpr-1.875d0*al1*                            &
     & (1.d0+al1*(al1*0.1875d0-0.5d0))                                  
      if(kentr.eq.0) go to 56 
      if(lst.eq.0) go to 53 
      st=cu(3)*(cu(4)+al1*(cu(27)+al1*(cu(5)*al1-cu(7)))) 
      st=st-cu(28)*qa*(cu(17)+al1*(cu(5)+cu(25)*al1)) 
      spl=x7-cu(4) 
   53 if(t.lt.tpar) go to 101 
      x8=x4+x5*(cu(30)+cu(31)*al1) 
      se=se+x8 
      if(lst.eq.0) go to 10 
      spl=spl-cu(15)*x8 
      st=st+x4*(cu(15)*alf+cu(34)+al1*(cu(32)+cu(33)*al1)) 
  101 if(lst.eq.0) go to 10 
      st=st*rg/t 
      spl=spl*rg 
   10 se=se*rg 
   56 pe=pi*pe 
      ee=eg*ee 
      if(lst.eq.0) go to 50 
      pt=rg*pl*pt 
      et=eg*et/t 
      ppl=ei*ppl 
      epl=eg*epl 
   50 go to 135 
!                                                                       
! *** addition the ion and black-body radiation components to eos       
   57 continue 
! ********************************************************************* 
      x=pi/zs 
      x1=eg/zs 
      fac=t*t 
      v=7.56471d-3*fac*fac 
!     v=7.56471d-3*t**4                                                 
      x3=v/cu(17) 
      p=x+pe+x3 
      beta=x3/p 
      x3=v/pl 
      e=x1+x3+ee 
      if(kentr.eq.0) go to 7 
      x6=cu(2)/as 
      x4=pl*(cu(35)/(t*sqrt(t))) 
      x4=log(x4) 
      x4=x6*(cu(24)-x4+scn) 
      x5=cu(36)*x3/t 
      s=x5+x4+se 
      sek=se/cu(2) 
      sk=s/cu(2) 
      if(lst.eq.0) go to 9 
      st=st+(cu(17)*x5+x1/t)/t 
      spl=spl-x6-x5 
      go to 45 
    7 if(lst.eq.0) go to 9 
   45 pt=pt+(x+cu(36)*v)/t 
      ppl=ppl+x/pl 
      et=et+(x1+cu(37)*x3)/t 
      epl=epl-x3 
! ********************************************************************* 
      x4=pt/(pl*et) 
      x5=pl/p 
      gam=x5*(ppl+t*x4*pt/pl) 
      da=x4/gam 
      cp=gam*et*(p/ppl) 
      dpe=x5*(epl+t*pt/pl)-1.d0 
      if(kentr.ne.0) then 
      dse=t*st/et-1.d0 
      dsp=-spl*(pl/pt)-1.d0 
      endif 
!                                                                       
!                 *** exit from epeos ***                               
    9 return 
!                                                                       
! *** further search for working region                                 
   51 if(plm.lt.y) go to 52 
      kzz=1 
  121 z=t/ep2 
      if(z.gt.0.57142857d0) go to 64 
      x=xep2*sq 
      if(kzz.eq.3) xz=cu(3)*x 
      go to 125 
   64 if(z.gt.1.07142857d0) go to 54 
      xz=7.306392d-2*z 
      x=xz+3.414385d-2 
      go to 125 
   54 if(z.gt.1.d1) go to 58 
      x=((4.77856d-2*z-1.41839d-2)*z+7.2001d-2)*z-7.20897d-3 
      if(kzz.eq.3) xz=((1.433568d-1*z-2.83678d-2)*z+7.2001d-2)*z 
      go to 125 
   58 fac=z*z*z 
      x=4.708d-2*fac 
!  58 x=4.708d-2*z**3                                                   
      if(kzz.eq.3) xz=cu(17)*x 
  125 go to (55,115,114),kzz 
   55 if(plm-x) 59,59,61 
!                                                                       
! *** expansion over half-integer fermi--dirac functions (nz=2)         
    2 sq=t*sqrt(t) 
   59 nzr=2 
      kenf=1 
      zp=plm/(cu(38)*sq) 
      x=psi 
      nit=0 
      kf=1 
      x1=cu(9)*al1 
      al2=al1*al1 
      al3=0.21875d0*al2 
      go to 26 
   11 z=zp-f12-x1*f32-al3*f52 
      v=f12s+x1*f32s+al3*f52s 
      dl=z/v 
      z1=x 
      if(z1.lt.0.d0)z1=-x 
      v=1.d0 
      if(dl.lt.0.d0) v=-1.d0 
      z=v*dl 
      if(z1.lt.0.1d0) go to 41 
      z=z/z1 
      if(z.gt.0.3d0) dl=0.3d0*v*z1 
   41 x=x+dl 
      nit=nit+1 
      if(z.lt.eit) go to 73 
      if(nit.lt.nitm) go to 26 
   66  write(6,5072) nzr 
       write(6,5073) dl,x,eit,nitm,jkk 
!..       print 72,nzr                                                  
!..       print 4072,dl,x,eit,nitm,jkk                                  
       stop 
   72  format('  iterations in region nzr=',i1,                         &
     & 1x,'do not converge!')                                           
 4072  format('  dx=',1p,d10.3,' x=',d10.3,' eit=',d10.3,               &
     & 1x,'nitm=',i3,' jkk=',i4,' *stop in epeos*')                     
 5072  format(10x,'epeos  ***  runaway of iterations in region nzr=',   &
     & i1,' *** stop*')                                                 
 5073  format(10x,1p,'dx=',d10.3,' x=',d10.3,' eit=',d10.3,' nitm=',i3, &
     & 1x,'jkk=',i4)                                                    
!                                                                       
   73 kf=2 
      go to 26 
  300 go to(31,9),kenf 
   31 psi=x 
      v=sq*al1 
      x=cu(19)*al1 
      z=cu(39)*v 
      z1=cu(3)*z/pl 
      z3=x1*f32 
      z4=x*f52 
      al4=9.375d-2*al2 
      y=al4*f72 
      pe=z*(f32+z4+y) 
      z5=x1*f52 
      y1=al3*f72 
      ee=z1*(f32+z5+y1) 
      if(kentr.eq.0) go to 34 
      z2=cu(41)*sq/pl 
      dl=cu(40)*psi 
      z10=cu(44)*psi 
      z11=cu(45)*al1 
      z6=z11*(f52-dl*f32) 
      al5=0.16875d0*al2 
      y3=0.77777777778d0*psi 
      y2=al5*(f72-y3*f52) 
      se=z2*(f32-z10*f12+z6+y2) 
   34 if(lst.eq.0) go to 35 
      pal=f12s+x1*f32s+al3*f52s 
      pap=zp/pal 
      pal=(cu(3)*zp+z3+cu(15)*al3*f52)/pal 
      z8=f32s+x1*f52s+al3*f72s 
      z7=z5+cu(15)*y1-pal*z8 
      et=(cu(27)*ee+z1*z7)/t 
      epl=z1*z8*pap-ee 
      z9=f32s+x*f52s+al4*f72s 
      pt=(cu(27)*pe+z*(z4+cu(15)*y-pal*z9))/t 
      ppl=z*z9*pap/pl 
      if(kentr.eq.0) go to 35 
      z9=f32s-z10*f12s-cu(44)*f12+z11*(f52s-dl*f32s-cu(40)*f32) 
      z9=z9+al5*(f72s-0.77777777778d0*f52-y3*f52s) 
      st=(cu(3)*se+z2*(z6+cu(15)*y2-pal*z9))/t 
      spl=z2*pap*z9-se 
   35 go to 135 
!                                                                       
! *** procedure of calculation of half-integer fermi--dirac functions   
      entry fd12f 
      kenf=2 
      kf=2 
      x=psi 
   26 if(x.ge.-1.d0) go to 21 
      z=exp(x) 
      v=ck(1)*z 
      f12=v*(1.d0+z*(a(1)+z*(a(2)+z*(a(3)+z*a(4))))) 
      f12s=v*(1.d0+z*(c(1)+z*(c(2)+z*(c(3)+z*c(4))))) 
      f32=ck(2)*z*(1.d0+z*(a(5)+z*(a(6)+z*(a(7)+z*a(8))))) 
      f52=ck(3)*z*(1.d0+z*(a(9)+z*(a(10)+z*(a(11)+z*a(12))))) 
      go to(30,12),kf 
   12 f72=a(13)*z*(1.d0+z*(a(14)+z*(a(15)+z*(a(16)+z*a(17))))) 
      go to 30 
   21 if(x.ge.-.1d0) go to 22 
      n=1 
      go to 14 
   22 if(x.ge.1.d0) go to 23 
      n=2 
      go to 14 
   23 if(x.ge.2.5d0) go to 24 
      n=3 
      go to 14 
   24 if(x.ge.4.5d0) go to 25 
      n=4 
   14 f12=a1(n)+x*(a2(n)+x*(a3(n)+x*a4(n))) 
      f12s=a2(n)+x*(df1(n)+x*df2(n)) 
      f32=b1(n)+x*(b2(n)+x*(b3(n)+x*(b4(n)+x*b5(n)))) 
      f52=c1(n)+x*(c2(n)+x*(c3(n)+x*(c4(n)+x*(c5(n)+x*c6(n))))) 
      go to(30,13),kf 
   13 f72=d1(n)+x*(d2(n)+x*(d3(n)+x*(d4(n)+x*(d5(n)+x*(d6(n)+x*d(n)))))) 
   30 f32s=1.5d0*f12 
      f52s=2.5d0*f32 
      if(kf.eq.1) go to 11 
      f72s=3.5d0*f52 
      go to 300 
   25 z=sqrt(x) 
      z1=z*x 
      z2=1.d0/(x*x) 
      f12=z1*(b(1)+(b(2)+(b(3)+b(4)*z2)*z2)*z2) 
      f12s=z*(c(5)+(c(6)+(c(7)+c(8)*z2)*z2)*z2) 
      z=z1*x 
      f32=z*(b(5)+(b(6)+(b(7)+b(8)*z2)*z2)*z2) 
      f52=x*z*(b(9)+(b(10)+(b(11)+b(12)*z2)*z2)*z2) 
      f32=f32+b(17) 
      f52=f52+b(18)+b(19)*x 
      go to(30,17),kf 
   17 f72=z*(b(13)+(b(14)+(b(15)+b(16)*z2)*z2)*z2)/z2 
      f72=f72+b(20)+x*(b(21)+b(22)*x) 
      go to 30 
!                                                                       
! *** search for working region is continued                            
   61 y=x*grcf 
      if(plm.lt.y) go to 62 
!                                                                       
! *** chandrasekhar's expansion for degenerate gas (nz=3)               
    3 nzr=3 
      x1=plm/cu(47) 
      if(psi.lt.3.d0) psi=3.d0 
      nit=0 
      x=psi*al1 
      x=sqrt(x*(x+cu(15))) 
      al2=al1*al1 
      z2=1.644934d0*al2 
      z4=1.894066d0*al2*al2 
   74 x2=x*x 
      x3=x2*x 
      x5=cu(17)*(z4/x3)/x2 
      x6=z2/x 
      x8=cu(15)*z2*x 
      dl=(x3*cu(49)+x8+x6+x5-x1)/(x3+x8-x6-cu(50)*x5) 
      x6=1.d0 
      if(dl.lt.0.d0) x6=-1.d0 
      z1=x6*dl 
      if(z1.gt.0.9d0) dl=0.9d0*x6 
      x=x*(cu(4)-dl) 
      nit=nit+1 
      if(z1.lt.eit) go to 71 
      if(nit.lt.nitm) go to 74 
      go to 66 
   71 x2=x*x 
      x4=cu(4)+x2 
      z=sqrt(x4) 
      y1=cu(4)+z 
         z1=x2/y1 
      psi=alf*z1 
      z3=x*z 
      x3=x*x2 
      z5=cu(15)*x2 
      x7=z5+cu(4) 
      if(x.gt.0.1d0) go to 174 
      x5=x2*x3 
      f0=x5*(1.6d0-x2*(fgs(1)-x2*(fgs(2)-x2*(fgs(3)-x2*fgs(4)))))*cu(49) 
      g0=x5*(0.8d0-x2*(fgs(5)-x2*(fgs(6)-x2*(fgs(7)-x2*fgs(8))))) 
      go to 175 
  174 x6=log(x+z) 
      f0=z3*(z5-cu(17))*cu(49)+x6 
      g0=x7*z3-x6-x3*cu(52) 
  175 x5=z4/x3 
      y5=z5-cu(4) 
      f2=z2*cu(51)*z3 
      f4=x5*cu(51)*y5*z 
      pe=cu(46)*(f0+f2+f4) 
      y2=cu(51)/y1 
      y4=z*y1 
      g2=z2*y2*x*(x7+y4) 
      g4=x5*cu(17)*y2*(cu(4)+y5*y4) 
      ee=cu(46)*(g0+g2+g4)/pl 
      if(kentr.eq.0) go to 75 
      y6=cu(48)/(t*pl) 
      se=y6*(f2+cu(15)*f4) 
   75 if(lst.eq.0) go to 76 
      z6=cu(54)*x5/x3 
      z7=x7-cu(15) 
      pap=x2+z2*z7/x2-z6 
      pal=cu(15)*(z2*x7+cu(55)*x5/x)/(x*pap) 
      pap=x1/pap 
      z9=cu(51)*x/z 
      z10=x1*z9 
      z11=cu(46)*pap/pl 
      ppl=z11*z10 
      y3=cu(46)/t 
      pt=y3*(cu(15)*(f2+cu(15)*f4)-z10*pal) 
      g0=z1*x2*cu(51) 
      v=cu(15)*(g2+cu(15)*g4) 
      g2=cu(51)*z2*((cu(4)+cu(17)*x7*y1)/y4-cu(15)) 
      g4=cu(53)*x5*(cu(37)-z)/(x*y4) 
      g0=g0+g2+g4 
      epl=z11*g0-ee 
      et=y3*(v-pal*g0)/pl 
      if(kentr.eq.0) go to 76 
      g4=(z2*x7+cu(55)*x5/x)*z9/x 
      spl=y6*pap*g4-se 
      st=(se+y6*(cu(37)*f4-pal*g4))/t 
   76 go to 135 
!                                                                       
! *** quadratures taken with the gauss method (nz=5)                    
    5 nzr=5 
      kkk=1 
         kk1=1 
      al3=alf*alf*alf 
         z11=plm*al3 
         x2=cu(15)*alf 
      z10=z11/cu(47) 
      kw=1 
         ku=10 
      go to 151 
  169 z=(g1m-z10)/g1mp 
      z3=1.d0 
      if(z.lt.0.d0) z3=-1.d0 
      z2=z*z3 
      z1=pc 
      if(pc.lt.0.d0) z1=-pc 
      if(z1.lt.1.d0) go to 181 
      z2=z2/z1 
      if(z2.gt.0.3d0) z=0.3d0*z1*z3 
  181 pc=pc-z 
      if(z2.gt.eit) go to 151 
      ku=15 
         kw=2 
      go to 151 
  170 z=1.44059d0/(al3*alf) 
      z1=z/pl 
        z2=z/cu(17) 
      pe=z2*gp 
         ee=z1*ge 
         hpr=cu(15)*g1/z10 
      if(kentr.eq.0) go to 182 
      z3=pc+alf 
         y3=g3+g31+cu(49)*gp-z3*g1m 
      se=z1*y3/t 
  182 if(lst.eq.0) go to 183 
      pap=g1m/g1mp 
         x=g1a1-g1a 
      pal=(cu(17)*g1m-alf*(x+cu(15)*g1p))/g1mp 
      x1=g2p1-g2p 
         y2=g2a+g2a1-cu(15)*g2p 
      ppl=ei*cu(49)*x1/g1mp 
      pt=pe*(cu(37)-(alf*y2+x1*pal)/gp)/t 
      y1=g4p1-g4p-x2*g1p 
      epl=ee*(pap*y1/ge-cu(4)) 
      et=ee*(cu(37)-(alf*(g4a1+g4a)+x2*(g1-g4p+alf*g1a-x2*g1p)+y1*pal)/ &
     &ge)/t                                                             
      if(kentr.eq.0) go to 183 
      y1=g3p1-g3p+cu(49)*x1-g1m-z3*g1mp 
      spl=se*(pap*y1/y3-cu(4)) 
      st=g3a1+g3a+y2*cu(49)-g1m-z3*(g1a1-g1a+cu(15)*g1p)-cu(15)*g3p 
      st=se*(cu(17)-(pal*y1+st*alf)/y3)/t 
  183 go to 135 
  151 kpg=0 
  152 wo=0.d0 
         w1=0.d0 
         w2=0.d0 
      wop=0.d0 
         w1p=0.d0 
        w2p=0.d0 
      woa=0.d0 
         w1a=0.d0 
        w2a=0.d0 
      if(pc.gt.pc2) go to 155 
      if(kkk.eq.0) go to 158 
      do 157 k=1,15 
  157 uac(k+65)=sqrt(uac(k)+x2) 
      kkk=0 
  158 if(pc.gt.pc1) go to 156 
      if(pc.lt.-4.4d1) go to 163 
      x=exp(pc) 
      do 161 k=1,ku 
  161 uac(k+80)=uac(k+15)/(uac(k+30)*x+1.d0) 
      do 162 k=1,5 
      z=wk1(k)*wk4(k) 
      wo=wo+z 
      wop=wop+z/(1.d0+aio(k)*x) 
      z=wk2(k)*wk5(k) 
      w1=w1+z 
      w1p=w1p+z/(1.d0+ai1(k)*x) 
      if(kw.eq.1) go to 162 
      z=wk6(k)*wk3(k) 
      w2=w2+z 
      if(lst.eq.0) go to 162 
      woa=woa+wk4(k)/wk1(k) 
      w1a=w1a+wk5(k)/wk2(k) 
      w2p=w2p+z/(1.d0+ai2(k)*x) 
      w2a=w2a+wk6(k)/wk3(k) 
  162 continue 
      wop=wop*x 
         w1p=w1p*x 
      wo=wo*x 
        w1=w1*x 
      if(kw.eq.1) go to 163 
      w2=w2*x 
      if(lst.eq.0) go to 163 
      w1a=w1a*x 
         woa=woa*x 
      w2a=w2a*x 
         w2p=w2p*x 
  163 g1=w1+alf*wo 
      g1p=w1p+alf*wop 
      if(kw.eq.1) go to 164 
      g2=w2+x2*w1 
         g3=w2+alf*(w1+g1) 
        g4=w2+alf*w1 
      if(lst.eq.0) go to 164 
      g1a=w1a+wo+alf*woa 
      g2p=w2p+x2*w1p 
         g3p=w2p+alf*(w1p+g1p) 
         g4p=w2p+alf*w1p 
      g2a=w2a+2.d0*w1+x2*w1a 
         g3a=w2a+w1+g1+alf*(w1a+g1a) 
      g4a=w2a+w1+alf*w1a 
  164 if(kpg.eq.1) go to 166 
      g11=g1 
         g1p1=g1p 
      if(kw.eq.1) go to 168 
      g21=g2 
         g31=g3 
        g41=g4 
      if(lst.eq.0) go to 168 
      g1a1=g1a 
         g2p1=g2p 
         g3p1=g3p 
         g4p1=g4p 
      g2a1=g2a 
         g3a1=g3a 
         g4a1=g4a 
  168 pc=-pc-x2 
         kpg=1 
      if(pc.gt.-4.4d1) go to 152 
         g1=0.d0 
         g2=0.d0 
         g3=0.d0 
         g4=0.d0 
         g1p=0.d0 
         g2p=0.d0 
         g3p=0.d0 
         g4p=0.d0 
         g1a=0.d0 
         g2a=0.d0 
         g3a=0.d0 
         g4a=0.d0 
  166 pc=-pc-x2 
      g1m=g11-g1 
      g1mp=g1p1+g1p 
      gp=g2+g21 
      ge=g4+g41+x2*g1 
  167 go to(169,170),kw 
  155 do 171 k=1,5 
      z4=xxi(k)-1.d0 
      z1=exp(pc*z4) 
      y1=pc*xxi(k) 
      z2=1.d0+z1 
      z3=x2+y1 
      y2=pc*aai(k)*sqrt(pc*z3)/z2 
      y4=cci(k)+pc 
      z=y4+x2 
      y6=bbi(k)*sqrt(y4*z) 
      wo=wo+y2+y6 
      if((lst.eq.0).and.(kw.ne.1)) go to 172 
      z5=1./pc 
         y3=0.5d0*xxi(k)/z3-z4*z1/z2+1.5d0*z5 
      z6=1.d0/y4 
         y5=0.5d0*(1.d0/z+z6) 
      wop=wop+y2*y3+y6*y5 
      z3=y2/z3 
         z=y6/z 
  172 y2=y2*y1 
         y6=y6*y4 
      w1=w1+y2+y6 
      y3=y3+z5 
        y5=y5+z6 
      w1p=w1p+y2*y3+y6*y5 
      if(kw.eq.1) go to 171 
      if(lst.eq.0) go to 173 
      woa=woa+z3+z 
      z3=z3*y1 
         z=z*y4 
      w1a=w1a+z3+z 
  173 y2=y2*y1 
         y6=y6*y4 
      w2=w2+y2+y6 
      if(lst.eq.0) go to 171 
      w2p=w2p+y2*(y3+z5)+y6*(y5+z6) 
      w2a=w2a+z3*y1+z*y4 
  171 continue 
      go to 163 
  156 if(kk1.eq.0) go to 191 
      do 197 k=1,5 
      wk6(k)=hhi(k) 
         wk7(k)=ggi(k) 
      wk4(k)=sqrt(vvi(k)+x2) 
  197 wk5(k)=sqrt(zzi(k)+x2) 
      do 195 k=1,24 
  195 abcg(k)=0.d0 
      k1=0 
      do 198 i=1,3 
      x=i 
         x1=x-0.5d0 
      do 190 k=1,5 
      k2=k1+k 
        k3=k2+65 
      z=(x+ggsi(k))/pc2 
         y=x1/zzi(k) 
      z1=wk7(k)*(z-cpp(2))*cpp(4) 
         z2=wk6(k)*(y-cpp(2))*cpp(4) 
      z4=xxi(k)*wk7(k)/wk4(k) 
        z5=wk6(k)/wk5(k) 
         z6=cpp(5)*(z4+z5) 
      asp(i)=asp(i)+abc(k2)*uac(k3)+z1*wk4(k)+z2*wk5(k)+z6*cpp(1) 
      z8=cpp(5)*(z4/(vvi(k)+x2)+z5/(zzi(k)+x2)) 
      aspa(i)=aspa(i)+abc(k2)/uac(k3)+z1/wk4(k)+z2/wk5(k)-z8*cpp(1) 
      k4=k2+15 
      z1=wk7(k)*(cpp(3)-z)*cpp(1) 
         z2=wk6(k)*(cpp(3)-y)*cpp(1) 
      bsp(i)=bsp(i)+abc(k4)*uac(k3)+z1*wk4(k)+z2*wk5(k)-z6 
      bspa(i)=bspa(i)+abc(k4)/uac(k3)+z1/wk4(k)+z2/wk5(k)+z8 
      k4=k4+15 
      csp(i)=csp(i)+abc(k4)*uac(k3) 
      cspa(i)=cspa(i)+abc(k4)/uac(k3) 
      k4=k4+15 
      gsp(i)=gsp(i)+abc(k4)*uac(k3) 
      gspa(i)=gspa(i)+abc(k4)/uac(k3) 
      wk6(k)=wk6(k)*zzi(k) 
  190 wk7(k)=wk7(k)*vvi(k) 
  198 k1=k1+5 
      kk1=0 
  191 z=pc-pc1 
         z1=2.d0*z 
         z2=1.5d0*z 
      wo=gsp(1)+z*(csp(1)+z*(bsp(1)+z*asp(1))) 
      w1=gsp(2)+z*(csp(2)+z*(bsp(2)+z*asp(2))) 
      wop=csp(1)+z1*(bsp(1)+z2*asp(1)) 
      w1p=csp(2)+z1*(bsp(2)+z2*asp(2)) 
      if(kw.eq.1) go to 163 
      w2=gsp(3)+z*(csp(3)+z*(bsp(3)+z*asp(3))) 
      if(lst.eq.0) go to 163 
      w2p=csp(3)+z1*(bsp(3)+z2*asp(3)) 
      woa=gspa(1)+z*(cspa(1)+z*(bspa(1)+z*aspa(1))) 
      w1a=gspa(2)+z*(cspa(2)+z*(bspa(2)+z*aspa(2))) 
      w2a=gspa(3)+z*(cspa(3)+z*(bspa(3)+z*aspa(3))) 
      go to 163 
!                                                                       
! *** relativistic asymptotics (nz=4)                                   
    4 nzr=4 
  520  ro=pl*cu(58) 
          hi=1.d0/emue 
          r1=ro*0.5d0*hi 
          pit=pi2*al1 
          pt2=pit*al1 
          pa=pt2-1.5d0 
      hu=psi*al1+1.d0 
      do 525 it=1,4 
          hu1=hu 
          hu2=hu*hu 
          hu =2.d0*(hu2*hu+r1)/(3.d0*hu2+pa) 
      if(abs(hu1/hu-1.d0).le.eit) go to 527 
  525  continue 
          fac=pa*.3333333333d0 
          r=r1*r1+fac*fac*fac 
!         r=r1**2+(pa*.3333333333d0)**3                                 
      x=(r1+sqrt(r))**(1.d0/3.d0) 
          hu=x-pa/(3.d0*x) 
  527  continue 
          hu2=hu*hu 
          psi=-7.77d2 
          if(al1.gt.1.d-8) psi=(hu-1.d0)/al1 
          pe=.25d0*(hu2*hu2+2.d0*pa*hu2+.46666666667d0*pt2*pt2-pt2) 
          ee=(3.d0*pe+0.5d0*(3.d0*hu2+pt2))/ro-hi 
          pe=pe+1.1550858d0 
          ee=ee-1.1550858d0/ro 
          if(lst.eq.0) go to 555 
          r=1.d0/(3.d0*hu2+pa) 
          r1=hu2+pa 
          pt=pit*(hu2-0.5d0+0.46666666667d0*pt2)-2.d0*pit*hu2*r1*r 
          r2=hi*hu*r 
          ppl=r2*r1 
          et=pit*(hu2*(1.d0-4.d0*pt2*r)+1.4d0*pt2-0.5d0)/ro 
          epl=(3.d0*(ppl+r2)-ee-hi) 
  555  ee=ee*cu(59) 
       pe=pe*cu(60) 
      if(kpar.ne.1) go to 558 
      hpr=0.d0 
      if(t.lt.5.96d0) go to 558 
      eta=alf*hu 
      if(eta.gt.6.d1) go to 558 
      hu1=exp(-eta) 
      al2=al1*al1 
      hpr=hu1*(1.2d1*al2-3.d0+hu1*((0.444444444d0*al2-1.d0)             &
     & *hu1-1.5d0*(al2-1.d0)))/(eta*r1)                                 
  558 continue 
          if(lst.eq.0)go to 556 
          pt=pt*cu(61) 
          ppl=ppl*cu(59) 
          epl=epl*cu(59) 
          et=et*cu(2) 
  556  if(kentr.eq.0) go to 557 
          y=cu(62) 
      se=y*al1*(hu2+.466666666667d0*pt2-.5d0)/pl 
          if(lst.eq.0) go to 557 
      spl=-se+2.d0*pi2*cu(2)*al1*r2 
          st=se/t+2.d0*y*pt2*(0.46666666667d0-2.d0*hu2*r)/(pl*cu(1)) 
  557  go to 135 
!                                                                       
! *** interpolation between perfect gas and expansion over              
! *** half-integer f-d functions                                        
   52 pl1=x*emue 
         pl2=y*emue 
        nzp1=1 
         nzp2=2 
!                                                                       
! *** interpolation over density                                        
   83 kin=ki 
      lst1=lst 
         lst=1 
         kk=0 
      dni=pl 
      if(lst1.eq.0) go to 81 
      kk=1 
      tni=t 
         t=t*(cu(4)+dst) 
   81 pl=pl1 
         nz=nzp1 
        ki=2 
      go to 90 
   77 pn1=pe 
         en1=ee 
        sn1=se 
         psn1=psi 
         hprn1=hpr 
      pnp1=ppl 
         enp1=epl/pl 
         snp1=spl/pl 
      pl=pl2 
         nz=nzp2 
        ki=3 
         nzr1=nzr 
      go to 90 
   78 if(kk.eq.2) go to 92 
      wv=pl2-pl1 
         wv1=cu(4)/wv 
      wv4=cu(15)*wv1 
         wv3=dni-pl1 
         wv5=wv1*wv3 
   92 x=epl/pl 
      x1=pe-pn1 
         x2=x1*wv4 
         x3=x1*wv1 
      z1=pnp1+ppl-x2 
         z2=x3-z1-pnp1 
         z1=wv1*z1 
      pe=pn1+wv3*(pnp1+wv5*(z2+wv3*z1)) 
      x1=ee-en1 
         x2=x1*wv4 
         x3=x1*wv1 
      y1=enp1+x-x2 
         y2=x3-y1-enp1 
         y1=wv1*y1 
      ee=en1+wv3*(enp1+wv5*(y2+wv3*y1)) 
      if(kentr.eq.0) go to 91 
      x1=se-sn1 
         x2=x1*wv4 
         x3=x1*wv1 
      v1=snp1+spl/pl-x2 
         v2=x3-v1-snp1 
         v1=wv1*v1 
      se=sn1+wv3*(snp1+wv5*(v2+wv3*v1)) 
   91 if(kk.eq.1) go to 82 
      if(kk.eq.2) go to 124 
  127 x1=(dni-pl1)*wv1 
      psi=(psi-psn1)*x1+psn1 
      hpr=(hpr-hprn1)*x1+hprn1 
      nzr=10*nzr1+nzr 
       pl=dni 
        plm=pl/emue 
      ki=kin 
         lst=lst1 
  134 nz=0 
         pi=ei*pl 
  135  go to(57,77,78,79,80,577,578), ki 
   82 pn2=pe 
         en2=ee 
        sn2=se 
      kk=2 
         t=tni 
      go to 81 
  124 x1=t*dst 
      pt=(pn2-pe)/x1 
         et=(en2-ee)/x1 
      if(kentr.eq.0) go to 126 
      st=(sn2-se)/x1 
      spl=snp1+wv5*(cu(15)*v2+cu(17)*wv3*v1) 
      spl=dni*spl 
  126 ppl=pnp1+wv5*(cu(15)*z2+cu(17)*wv3*z1) 
      epl=enp1+wv5*(cu(15)*y2+cu(17)*wv3*y1) 
      epl=dni*epl 
      go to 127 
!                                                                       
! *** interpolation between degenerate gas and expansion over           
! *** half-integer fermi-dirac functions                                
   62 pl1=x*emue 
         pl2=y*emue 
        nzp1=2 
         nzp2=3 
      go to 83 
!                                                                       
! *** interpolation over temperature                                    
  128 kit=ki 
      lst2=lst 
         lst=1 
         kkt=0 
      dnt=t 
      if(lst2.eq.0) go to 129 
      kkt=1 
         plni=pl 
        pl=pl*(cu(4)+dst) 
  129 t=t1 
         nz=nz1 
         ki=4 
      go to 90 
   79 pnt=pe 
         ent=ee 
         snt=se 
         psnt=psi 
         hprnt=hpr 
         pntt=pt 
         entt=et 
         sntt=st 
         t=t2 
         nz=nz2 
         ki=5 
         nzrt=nzr 
      go to 90 
   80 if(kkt.eq.2) go to 93 
      vw=t2-t1 
         vw1=cu(4)/vw 
      vw4=cu(15)*vw1 
         vw3=dnt-t1 
        vw5=vw1*vw3 
   93 x1=pe-pnt 
         x2=x1*vw4 
         x3=x1*vw1 
      z1=pntt+pt-x2 
         z2=x3-z1-pntt 
         z1=vw1*z1 
      pe=pnt+vw3*(pntt+vw5*(z2+vw3*z1)) 
      x1=ee-ent 
         x2=x1*vw4 
         x3=x1*vw1 
      y1=entt+et-x2 
         y2=x3-y1-entt 
         y1=vw1*y1 
      ee=ent+vw3*(entt+vw5*(y2+vw3*y1)) 
      if(kentr.eq.0) go to 94 
      x1=se-snt 
         x2=x1*vw4 
         x3=x1*vw1 
      v1=sntt+st-x2 
         v2=x3-v1-sntt 
         v1=vw1*v1 
      se=snt+vw3*(sntt+vw5*(v2+vw3*v1)) 
   94 if(kkt.eq.1) go to 130 
      if(kkt.eq.2) go to 131 
  133 x1=(dnt-t1)*vw1 
      psi=(psi-psnt)*x1+psnt 
      hpr=(hpr-hprnt)*x1+hprnt 
      nzr=10*nzrt+nzr 
      t=dnt 
         alf=cu(1)/t 
         al1=cu(4)/alf 
         ei=t*rg 
      eg=cu(3)*ei 
         ki=kit 
         lst=lst2 
      go to 134 
  130 pnt2=pe 
        ent2=ee 
         snt2=se 
      kkt=2 
         pl=plni 
      go to 129 
  131 x1=cu(4)/dst 
      ppl=x1*(pnt2-pe)/pl 
         epl=x1*(ent2-ee) 
      if(kentr.eq.0) go to 132 
      spl=(snt2-se)*x1 
      st=sntt+vw5*(cu(15)*v2+cu(17)*vw3*v1) 
  132 pt=pntt+vw5*(cu(15)*z2+cu(17)*vw3*z1) 
      et=entt+vw5*(cu(15)*y2+cu(17)*vw3*y1) 
      go to 133 
      END                                           
                                                                        
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         INVEOS                                               
!    MODULE:       INVEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         5/23/90                                              
!                  Bug fixed on (5/24/90)                               
!                                                                       
!                                                                       
!    CALL LINE:    CALL INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,XPREV
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  T_OLD = INITIAL GUESS AT THE TEMPERATURE             
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                  IFLAG = 1 --> INPVAR IS TEMPERATURE                  
!                          2 --> INPVAR IS INTERNAL ENERGY              
!                          3 --> INPVAR IS ENTROPY (NOT IMPLEM)         
!                                                                       
!    OUTPUTS       EOSFLG = 1 --> "NO NUCLEI" EOS                       
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS FOR DENSITIES ABOVE NUCLEAR  
!                  XPREV = PREVIOUS VALUE OF X                          
!                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY            
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M1D.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,            &
     &                  FORFLG,SF,XPREV,P_PREV)                         
!                                                                       
!                                                                       
!                                                                       
                                                                        
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
!                         Local variables                               
      INCLUDE 'eos_m4a.inc' 
!                                                                       
      DOUBLE PRECISION INP_V, INP_VO, INP_VN, UFTN, DUFTN, DT 
      DOUBLE PRECISION T_OLD, T_NEW, T_TEMP, T_LB, T_UB, PERDIF 
      DOUBLE PRECISION eps, prec, var, t_old0, al, fa, bl, fb, cl,      &
     &                 fc, d, e, tm, tol1, s, p, ql, r                  
      INTEGER LOOP, SF, NEW_F 
      INTEGER itmax, IFLAG0 
!                                                                       
      data itmax/300/ , eps/1.d-20/ 
      data prec/1.d-14/ 
!                                                                       
      RSFLAG = 1 
!                         Input is the temperature; call the EOS        
!                         normally and then return                      
      IF(IFLAG.EQ.1) THEN 
        CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,               &
     & XPREV,P_PREV)                                                    
        T_OLD = INPVAR(1) 
        RETURN 
      ENDIF 
!                                                                       
!--the calling variable is either internal energy or entropy            
!                                                                       
      IFLAG0=IFLAG 
      var=INPVAR(1) 
      t_old0=T_OLD 
!                                                                       
!--find lower bound in T                                                
!                                                                       
      INPVAR(1)=0.5*t_old0 
      al=INPVAR(1) 
      IFLAG=1 
      CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,                 &
     &             XPREV,P_PREV)                                        
      IF(IFLAG0.eq.2)THEN 
         fa=UTOT - var 
      ELSE 
         fa=STOT - var 
      END IF 
!                                                                       
!--find upper bound in T                                                
!                                                                       
      INPVAR(1)=2.*t_old0 
      bl=INPVAR(1) 
      IFLAG=1 
      CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,                 &
     &             XPREV,P_PREV)                                        
      IF(IFLAG0.eq.2)THEN 
         fb=UTOT - var 
      ELSE 
         fb=STOT - var 
      END IF 
!                                                                       
!--root has been bracketed, iterate to find solution                    
!                                                                       
   10 continue 
      fc=fb 
      do i=1,itmax 
!                                                                       
!--rename a,b,c and adjust bounding interval                            
!                                                                       
         if(fb*fc.gt.0.) then 
            cl=al 
            fc=fa 
            d=bl-al 
            e=d 
         end if 
         if(dabs(fc).lt.dabs(fb)) then 
            al=bl 
            bl=cl 
            cl=al 
            fa=fb 
            fb=fc 
            fc=fa 
         end if 
!                                                                       
!--check for convergence                                                
!                                                                       
         tm=0.5*(cl-bl) 
         tol1=2.*eps*dabs(bl) 
         if(dabs(tm).lt.tol1.or.dabs(fb/var).lt.prec) then 
            T_OLD=bl 
            return 
         end if 
!                                                                       
!--attempt inverse quadratic interpolation                              
!                                                                       
         if(dabs(e).ge.tol1.and.dabs(fa).gt.dabs(fb)) then 
            s=fb/fa 
            if(al.eq.cl) then 
               p=2.*tm*s 
               ql=1.-s 
            else 
               ql=fa/fc 
               r=fb/fc 
               p=s*(2.*tm*ql*(ql-r)-(bl-al)*(r-1.)) 
               ql=(ql-1.)*(r-1.)*(s-1.) 
            end if 
!                                                                       
!--check wether in bound                                                
!                                                                       
            if(p.gt.0.) ql=-ql 
            p=dabs(p) 
!                                                                       
!--accept or refuse interpolation                                       
!                                                                       
            if(2.*p.lt.dmin1(3.*tm*ql-dabs(tol1*ql),dabs(e*ql))) then 
!                                                                       
!--accept interpolation                                                 
!                                                                       
               e=d 
               d=p/ql 
            else 
!                                                                       
!--interpolation failed use bisection                                   
!                                                                       
               d=tm 
               e=d 
            end if 
!                                                                       
!--bound decreasing to slowly use bisection                             
!                                                                       
         else 
            d=tm 
            e=d 
         end if 
!                                                                       
!--move last guess to a                                                 
!                                                                       
         al=bl 
         fa=fb 
!                                                                       
!--evalue new trial point                                               
!                                                                       
         if(dabs(d).gt.tol1) then 
            bl=bl + d 
         else 
            bl=bl + dsign(tol1,tm) 
         end if 
         INPVAR(1)=bl 
         IFLAG=1 
         CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,              &
     &                XPREV,P_PREV)                                     
         IF(IFLAG0.eq.2)THEN 
            fb=UTOT - var 
         ELSE 
            fb=STOT - var 
         END IF 
!                                                                       
      enddo 
!                                                                       
!--did not converge write(iprint,error message                          
!                                                                       
      write(*,*)'iterations did not converge!' 
!                                                                       
      return 
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EOS4B.FOR                                            
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EOS_M4B                                              
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/15/90 Modified from model 4A to include the       
!                  phase boundary cutoffs and Maxwell construction      
!                  boundaries.                                          
!                  7/13/90 Modified from model 1-d to include Maxwell   
!                  construction                                         
!                  5/25/90  MODEL 1D   (RELEASE # 1.2)                  
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL EOS_M4A(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,    
!                  XPREV,P_PREV)                                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                  IFLAG = 1 --> INPVAR IS TEMPERATURE                  
!                          2 --> INPVAR IS INTERNAL ENERGY (NOT IMPLEM) 
!                          3 --> INPVAR IS ENTROPY (NOT IMPLEMENTED)    
!                  FFLAG = "FORCING FLAG"  0 --> NO FORCING             
!                                          1 --> FORCE A PARTICULAR     
!                                                SCHEME TO BE USED      
!                                                                       
!                                                                       
!    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B           
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS (includes alpha's)           
!                  XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON     
!                          FIRST CALL)                                  
!                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY (MUST BE   
!                          SUPPLIED ON FIRST CALL)                      
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EOS_M4A(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,SSFLAG,    &
     &                   XPREV,P_PREV)                                  
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION OUTVAR(4) 
!                                                                       
!                                                                       
!                       This include file contains all variable         
!                       declarations.  NOTE:: no implicit typing        
!                       scheme is followed in this code; if you         
!                       have any doubt as to a variables type CHECK     
!                       IT!!!!.  Also note that ALL variables are       
!                       declared explicitly.                            
!                                                                       
      INCLUDE 'eos_m4a.inc' 
!                                                                       
!                                                                       
!                         Set the "switch" flag to zero                 
      SWTFLG = 0 
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options should go through 
!                         INVEOS untill further notice)                 
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                         If the "forcing" flag is set then skip        
!                         the EOS determination logic and go straight   
!                         to the EOS determined by EOSFLG               
      IF(FFLAG.EQ.1) THEN 
        GOTO 10 
      ELSE 
!                         Otherwise let the EOS logic module determine  
!                         the correct EOS to use                        
        CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG) 
      ENDIF 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                        Try NUCEOS first and if not successfull        
!                        then try bulk EOS                              
   10 CONTINUE 
      IF(EOSFLG.EQ.1) THEN 
!                                                                       
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR) 
          OUTVAR(1) = INPVAR(1) 
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
          IF(SSFLAG.NE.1) THEN 
            OUTVAR(2) = 0.155 
            OUTVAR(3) = -15.0 
            OUTVAR(4) = -20.0 
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
          ENDIF 
!                                                                       
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN 
!                    Set EOS flag to full scheme                        
          EOSFLG = 2 
!                                                                       
!                    Else if fraction of nuclei is less than the minimum
!                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE 
          IF(FFLAG.NE.1) THEN 
!                                                                       
            CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
            IF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN 
              EOSFLG = 1 
              WRITE(*,*) 'A2 failed at try = ',T,BRYDNS,YE 
              GOTO 999 
            ENDIF 
!                                                                       
!                    Set nuclei to bulk EOS                             
            EOSFLG = 3 
!                    Save value of proton fraction                      
            P_PREV = YE*BRYDNS 
!                                                                       
            GOTO 999 
!                                                                       
          ELSE 
            IF(NF_FLG.EQ.1)                                             &
     &          WRITE(*,*) 'NUC failed at t,rho = ',t,brydns            
            GOTO 999 
          ENDIF 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!                                                                       
!          End of NUCEOS--BULK EOS calculations                         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                            CALCULATE FULL EOS (INCLUDING NUCLEI)      
      IF(EOSFLG.EQ.2) THEN 
!                                                                       
!                    Call the nuclear EOS                               
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
!ccc          WRITE(*,*) ' EOS_M4A:: r.i.gs.'                           
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR) 
          OUTVAR(1) = INPVAR(1) 
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
          IF(SSFLAG.NE.1) THEN 
            OUTVAR(2) = 0.155 
            OUTVAR(3) = -15.0 
            OUTVAR(4) = -20.0 
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
          ENDIF 
!                                                                       
!                                                                       
!                                                                       
          IF(SSFLAG.NE.1) THEN 
!ccc            WRITE(*,*) '     r.i.gs. failure @ try: ',inpvar        
            GOTO 999 
          ELSE 
            INPVAR(2) = OUTVAR(2) 
            INPVAR(3) = OUTVAR(3) 
            INPVAR(4) = OUTVAR(4) 
          ENDIF 
!                    Otherwise quit and return                          
        ELSEIF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN 
          GOTO 999 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    If fraction of heavies is greater than the minimum 
!                    parameter, then this EOS is OK                     
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN 
!                    Set EOS flag to full scheme                        
          EOSFLG = 2 
!                                                                       
!                    Else if fraction of nuclei is less than the minimum
!                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE 
!                    If the forcing flag is not set                     
          IF(FFLAG.NE.1) THEN 
!                    Set nuclei to no nuclei EOS                        
            EOSFLG = 3 
!                    Set flag to indicate switch is being made          
            SWTFLG = 1 
!                                                                       
            WRITE(*,*) ' NUCEOS failed at try =',t,brydns,ye 
            WRITE(*,*) ' where it shouldnt have; Bulk EOS was used' 
            WRITE(*,*) ' IV = ',INPVAR 
            WRITE(*,*) ' ' 
!                                                                       
!                    Branch to bulk EOS                                 
            GOTO 50 
!                                                                       
!                    Otherwise since forcing flag is set then declare   
!                    a failure and return                               
          ELSE 
!                      If the failure message flag is set then announce 
!                      the failure                                      
            IF(NF_FLG.EQ.1)                                             &
     &          WRITE(*,*) 'NUC failed at t,r = ',t,brydns              
            GOTO 999 
          ENDIF 
        ENDIF 
!                                                                       
      ENDIF 
!                              END OF FULL EOS CALULATIONS              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                              CALCULATE BULK EOS                       
   50 CONTINUE 
      IF(EOSFLG.EQ.3) THEN 
!                                                                       
        CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
        IF((SSFLAG.EQ.0).AND.(FFLAG.EQ.1).AND.(NF_FLG.EQ.1)) THEN 
          WRITE(*,*) 'A1 failed at t,rho = ',t,brydns 
          GOTO 999 
        ENDIF 
!                           If this EOS was used as a result of the     
!                           nuclear EOS failing then set the            
!                           success flag to indicate a warning          
        IF(SWTFLG.EQ.1) THEN 
          SSFLAG = 2 
        ENDIF 
!                                                                       
!                           Save the value of the proton fraction       
        P_PREV = YE*BRYDNS 
!                                                                       
        GOTO 999 
!                                                                       
      ENDIF 
!                END OF BULK EOS CALCULATIONS                           
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                              CALCULATE VIA MAXWELL CONSTRUCTION       
      IF(EOSFLG.EQ.4) THEN 
!                                                                       
        CALL MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                 Save the value of the proton fraction                 
        P_PREV = YE*BRYDNS 
!                                                                       
!                 If Maxwell EOS failed then announce the failure       
        IF(SSFLAG.NE.1) THEN 
          WRITE(*,*) ' MAXWEL failed at try = ' 
          WRITE(*,*) T,BRYDNS,YE 
        ENDIF 
!                                                                       
          GOTO 999 
!                                                                       
      ENDIF 
!                END OF MAXWELL CONSTRUCTION CALCULATIONS               
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         NUCEOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       NUCEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/13/90 Modified from model 1-d                      
!                                                                       
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!    CALL LINE:    CALL NUCEOS(INPVAR,YE,BRYDNS,X_PREV,SSFLAG)          
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:      XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON     
!                  FIRST CALL)                                          
!                  SSFLAG = SUCCESS FLAG 0 --> FAILURE                  
!                                        1 --> SUCCESS                  
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                       Function type declarations                      
!                                                                       
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFI, FHALFO 
      double precision fhalf 
!                                                                       
      DOUBLE PRECISION ZNG, ZPG 
      INTEGER TCFLAG, ftflag 
!                                                                       
      INTEGER KKI,LLI 
      DOUBLE PRECISION RESULT(5), R_CHECK(5) 
      double precision a_tmp(5,5) 
      DOUBLE PRECISION NI_MIN 
!                                                                       
      integer cflag, schflg 
      double precision dtst1, dtst2 
      double precision break, dnsi, dtmp8 
      double precision dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6,dtmp7 
!c      double precision tbsph, tbph, tbnh, tbspl, tbpl, tbnl           
!c      double precision dbspdx, dbpdx, dbndx, dbspdu, dbpdu, dbndu     
!c      double precision tsgl, tsgh, thl, thh, dsgdx, dhfdx, ds2dx,dzdx 
!c      double precision dpt1dx, dpt2dx                                 
!                                                                       
!                                                                       
!                         Set the scheme flag to zero                   
      SCHFLG = 0 
!                                                                       
!                                                                       
    5 CONTINUE 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options are not implemente
!                         in this version)                              
      T = INPVAR(1) 
      NSUBI = INPVAR(2) 
      ETA_PO = INPVAR(3) 
      ETA_NO = INPVAR(4) 
!                                                                       
!                                                                       
!                         Calc the quantum concentration of nucleons    
      NQ = 2.36D-4*T**1.5 
!                                                                       
!                         Calc the Fermi integral coefficent            
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
      KQ = ((T/UQ)**2.5)/(2.0*PI**2) 
!                                                                       
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2)) 
!                                                                       
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ) 
!                                                                       
      IF(ETA_PO.GE.ETAMAX) ETA_PO = ETAMAX-0.1 
      IF(ETA_NO.GE.ETAMAX) ETA_NO = ETAMAX-0.1 
      NI_MIN = DMAX1(4.5D-2,BRYDNS) 
      IF(NSUBI.LT.NI_MIN) NSUBI = NI_MIN+1.0D-3 
!                                                                       
      TCFLAG = 0 
!                                                                       
      cflag = 0 
!                                                                       
      NEWFLG = 1 
!                                                                       
!                    Start Newton-Raphson iteration here                
!                                                                       
!                                                                       
      DO 30 I=1,MAXIT,1 
!                                                                       
        IT_NUM = I 
!                       Set the "Negative" flag                         
        NGFLAG = 0 
!                                                                       
!                                                                       
!                                                                       
        NNOUT = MQ*F_1_2(ETA_NO)/(2.0*PI**2) 
        NPOUT = MQ*F_1_2(ETA_PO)/(2.0*PI**2) 
!                                                                       
        NOUT = NPOUT+NNOUT 
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
!                                                                       
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
!                                                                       
        MUALFA = 2.0*MUN_O+2.0*MUP_O+BALPHA-BPROUT*V_ALFA 
!                                                                       
        IF(ABS(MUALFA/T).LT.30.0) THEN 
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
        ELSEIF((MUALFA/T).LT.-30.0) THEN 
          ALFDNS = 0.0 
        ELSE 
          ALFDNS = 8.0*NQ*DEXP(3.0D1) 
        ENDIF 
!                                                                       
!                                                                       
!                   These statements take out the alfas if the          
!                   alpha particle enable flag is not set               
        IF(ALFLAG.NE.1) THEN 
          ALFDNS = 0.0 
          MUALFA = -300.0 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                                                                       
        BPRALF = ALFDNS*T 
!                                                                       
!---------------------------------------------------                    
!                                                                       
!                                                                       
!             Calculate fraction of space occupied by nuclei            
        U_NUC = (BRYDNS-EXALFA*NOUT-4.0*ALFDNS)/                        &
     &        (NSUBI-EXALFA*NOUT-4.0*ALFDNS)                            
!                                                                       
!                                                                       
!            Is volume occupied by nuclei within acceptable limits?     
!c        IF((U_NUC.LT.0.0).OR.((U_NUC-1.0).GT.-1.0E-20)) THEN          
        IF((U_NUC.LT.0.0).OR.(U_NUC.GT.0.996)) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!            Volume exclusion factor due to nuclei                      
        EXCLU = 1.0-U_NUC 
!                                                                       
!                                                                       
!            If calculated nucleon and alfa densities are larger        
!            than the baryon density then reduce the eta's              
        IF((EXCLU*EXALFA*NOUT+EXCLU*4.0*ALFDNS).GT.BRYDNS) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!            Calculate the internal (inside nuclei) proton fraction     
!                                                                       
        X = (BRYDNS*YE-(1.0-U_NUC)*(EXALFA*NPOUT+2.0*ALFDNS))/          &
     &    (U_NUC*NSUBI)                                                 
        COMPX = 1.0-X 
!                                                                       
!                                                                       
!            Is X within reasonable (but not necessarily correct)       
!            limits? (YE may not be the lower bound on X !!!)           
!ccc        X_MIN = DMAX1(1.0D-2,(YE-0.05))                             
        X_MIN = DMAX1(1.0D-2,(0.8*YE)) 
!c        x_min = 0.95*ye                                               
        IF((X.LT.X_MIN).OR.(X.GT.0.6)) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                     Calculate critical temperature & its X derivative 
        TSC_12 = 87.76*((COMP/375.0)**0.5)*((0.155/NSUBS)**OVR3) 
!                                                                       
!ccdebug      tsc_12 = 1.0d8                                            
!                                                                       
        TSUBC = TSC_12*X*COMPX 
        DTCDX = TSC_12*(1.0-2.0*X) 
        DTCDXX = -2.0*TSC_12 
!                                                                       
!c        tsubc = tsc_12*0.25                                           
!c        dtcdx = 0.0                                                   
!c        dtcdxx = 0.0                                                  
!                                                                       
!                                                                       
!C        TSUBC = 80.0*X*COMPX                                          
!C        DTCDX = 80.0*(1.0-2.0*X)                                      
!C        DTCDXX = -160.0                                               
!                                                                       
!                     If the X is such that T is greater than the       
!                     critical temperature then fix NSUBI so that       
!                     it lies in the bounds of acceptable parameter     
!                     space                                             
        ftflag = 0 
        IF((T.GT.TSUBC).AND.(SCHFLG.EQ.0)) THEN 
!                       If this is an initial guess, then lower         
!                       NSUBI untill we get a good X                    
          IF(NEWFLG.EQ.1) THEN 
!c        write(*,*) ' nuc exceeded Tc'                                 
!c        write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc                   
            ZNG = 2.0*(PI**2)*BRYDNS*(1.0-0.1*YE)/(1.5*MQ) 
            ZPG = 2.0*(PI**2)*BRYDNS*0.1*YE/(1.5*MQ) 
            IF(TCFLAG.NE.1) THEN 
              TCFLAG = 1 
              ETA_PO = FINV12(ZPG)-0.0 
              ETA_NO = FINV12(ZNG)-0.0 
            ELSE 
              ETA_PO = ETA_PO-2.0/T 
              ETA_NO = ETA_NO-2.0/T 
              NSUBI = DMAX1(0.9*NSUBI,5.1D-2) 
            ENDIF 
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '1',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
            GOTO 30 
          ELSE 
!                       Otherwise go back and cut the stepsize in       
!                       half since it was obviously too big             
            NGFLAG = 1 
            GOTO 29 
          ENDIF 
        ELSEIF((T.GT.TSUBC).AND.(SCHFLG.EQ.1)) THEN 
          ftflag = 1 
          tsubc = 80.0*(0.25+0.5*ye)*(0.75-0.25*ye) 
!                                                                       
        ENDIF 
!                                                                       
!                                                                       
        R_0 = (0.75/(PI*NSUBS))**OVR3 
        Q = (384.0*PI*(R_0**2)*SIG_S/SYM_S)-16.0 
!                                                                       
!                        Calculate surface functions of the internal    
!                        (nuclear) proton fraction, X                   
        SIGMA = 1.0/(Q+1.0/(X**3)+1.0/(COMPX**3)) 
        OVRX4 = (1.0/X**4)-(1.0/COMPX**4) 
        DSIGDX = 3.0*(SIGMA**2)*OVRX4 
        SIGSGP = DSIGDX/SIGMA 
        SIGSG2 = 18.0*(SIGMA**2)*OVRX4**2-12.0*SIGMA*((1.0/X**5)+       &
     &  (1.0/COMPX**5))                                                 
!                                                                       
!                        If T is less than critical temp then           
        IF(T.LT.TSUBC) THEN 
!                        Calculate the surface energy temperature factor
!                        and its X and T derivatives                    
!--this did not work because H could become 0 due to round-off          
!          H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4                        
!--this is better                                                       
          H = (1.0-(T/TSUBC)**2)**2 
          HPRIM = -4.0*T/(TSUBC**2)+4.0*((T/TSUBC)**3)/TSUBC 
          HPPRIM = -4.0/(TSUBC**2)+12.0*(T**2)/(TSUBC**4) 
          DHDX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDX 
          DHDXX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDXX+             &
     &    4.0*(-3.0*T**2/TSUBC**4+5.0*T**4/TSUBC**6)*(DTCDX**2)         
          HX = DHDX/H 
          DHDTDX = 8.0*(T/TSUBC**3-2.0*(T**3)/TSUBC**5)*DTCDX 
!                                                                       
!                                                                       
!                        X independent version of TZERO                 
!          TZERO = 0.25*TSC_12                                          
!          DTZDX = 0.0                                                  
!          DTZDXX = 0.0                                                 
!                        X dependent version of TZERO                   
!          TZERO = TSUBC                                                
!          DTZDX = DTCDX                                                
!          DTZDXX = DTCDXX                                              
!                                                                       
!                                                                       
!                                                                       
!                        Coulomb liquid correction factors and their    
!                        derivatives                                    
!          W = 1-(T/TZERO)**2                                           
!          DWDX = 2.0*(T**2)*DTZDX/(TZERO**3)                           
!          DWDT = -2.0*T/(TZERO**2)                                     
!          DWDTDX = 4.0*T*DTZDX/(TZERO**3)                              
!          DWDXDX = 2.0*(T**2)*                                         
!     1    (DTZDXX/(TZERO**3)-3.0*(DTZDX**2)/(TZERO**4))                
!          DWDTDT = -2.0/(TZERO**2)                                     
!                                                                       
          w = 1.0 
          dwdt = 0.0 
          dwdx = 0.0 
          dwdtdx = 0.0 
          dwdxdx = 0.0 
          dwdtdt = 0.0 
!                                                                       
!                                                                       
!                                                                       
!                        Calc lattice factor & derivatives & products   
!                                                                       
          EXCLU = 1.0-U_NUC 
          COMPU = 1.0-U_NUC 
!                                                                       
          DU = DMAX1(1.0D-15, (1.0-1.5*W*U_NUC**OVR3+0.5*U_NUC)) 
          DMU = DMAX1(1.0D-15,(1.0-1.5*W*(1.0-U_NUC+1.0E-20)**OVR3+     &
     &    0.5*(1.0-U_NUC)))                                             
!                                                                       
          DUP = -0.5*W*U_NUC**M2OVR3+0.5 
          DMUP =-0.5*W*(1.0-U_NUC+1.0E-20)**M2OVR3+0.5 
          DUPP = OVR3*W*((U_NUC+1.0D-20)**M5OVR3) 
          DMUPP = OVR3*W*((1.0-U_NUC)+1.0E-20)**M5OVR3 
!                                                                       
!                Derivatives w.r.t. T                                   
!                                                                       
          DUT = -1.5*DWDT*U_NUC**OVR3 
          DMUT = -1.5*DWDT*(1.0-U_NUC+1.0E-20)**OVR3 
          DUPT = -0.5*DWDT*U_NUC**M2OVR3 
          DMUPT = -0.5*DWDT*(1.0-U_NUC+1.0E-20)**M2OVR3 
!                                                                       
!                Derivatives w.r.t. X                                   
!                                                                       
          DUX = -1.5*DWDX*U_NUC**OVR3 
          DMUX = -1.5*DWDX*(1.0-U_NUC+1.0E-20)**OVR3 
          DUPX = -0.5*DWDX*U_NUC**M2OVR3 
          DMUPX = -0.5*DWDX*(1.0-U_NUC+1.0E-20)**M2OVR3 
!                                                                       
!                Second derivatives w.r.t. X                            
!                                                                       
          DUXX = -1.5*DWDXDX*U_NUC**OVR3 
          DMUXX = -1.5*DWDXDX*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                Second derivatives w.r.t. T                            
!                                                                       
          DUTT = -1.5*DWDTDT*U_NUC**OVR3 
          DMUTT = -1.5*DWDTDT*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                Second derivatives w.r.t. X & T                        
!                                                                       
          DUXT = -1.5*DWDTDX*U_NUC**OVR3 
          DMUXT = -1.5*DWDTDX*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                                                                       
          TMP1 = (U_NUC**2)+(COMPU**2)+0.6*(U_NUC*COMPU)**2 
          TMP1P = 4.0*U_NUC-2.0+                                        &
     &    2.0*0.6*(U_NUC*COMPU**2-COMPU*U_NUC**2)                       
          TMP1PP = 4.0+2.0*0.6*(COMPU**2-4.0*U_NUC*COMPU+U_NUC**2) 
!                                                                       
          TMP2 = COMPU*(DU**OVR3) 
          TMP2P = -1.0*DU**OVR3+OVR3*COMPU*(DU**M2OVR3)*DUP 
          TMP2PP = -OVR23*(DU**M2OVR3)*DUP-OVR29*COMPU*                 &
     &    (DU**M5OVR3)*DUP**2+OVR3*COMPU*(DU**M2OVR3)*DUPP              
!                                                                       
          TMP2T = OVR3*COMPU*(DU**M2OVR3)*DUT 
          TMP2X = OVR3*COMPU*(DU**M2OVR3)*DUX 
          TMP2XX = OVR3*COMPU*(DU**M2OVR3)*DUXX+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUX**2)                   
          TMP2TT = OVR3*COMPU*(DU**M2OVR3)*DUTT+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUT**2)                   
          TMP2XT = OVR3*COMPU*(DU**M2OVR3)*DUXT+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUX*DUT                    
          TMP2PT = -OVR3*(DU**M2OVR3)*DUT+                              &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUT+                   &
     &        OVR3*COMPU*(DU**M2OVR3)*DUPT                              
          TMP2PX = -OVR3*(DU**M2OVR3)*DUX+                              &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUX+                   &
     &        OVR3*COMPU*(DU**M2OVR3)*DUPX                              
!                                                                       
!                                                                       
!                                                                       
          TMP3 = U_NUC*(DMU**OVR3) 
          TMP3P = (DMU**OVR3)-OVR3*U_NUC*(DMU**M2OVR3)*DMUP 
          TMP3PP = -OVR23*(DMU**M2OVR3)*DMUP-OVR29*U_NUC*               &
     &    (DMU**M5OVR3)*(DMUP**2)+OVR3*U_NUC*(DMU**M2OVR3)*DMUPP        
!                                                                       
          TMP3T = OVR3*U_NUC*(DMU**M2OVR3)*DMUT 
          TMP3X = OVR3*U_NUC*(DMU**M2OVR3)*DMUX 
          TMP3XX = OVR3*U_NUC*(DMU**M2OVR3)*DMUXX+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUX**2)                 
          TMP3TT = OVR3*U_NUC*(DMU**M2OVR3)*DMUTT+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUT**2)                 
          TMP3XT = OVR3*U_NUC*(DMU**M2OVR3)*DMUXT+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*DMUX*DMUT                 
          TMP3PT = OVR3*(DMU**M2OVR3)*DMUT-OVR3*M2OVR3*U_NUC*           &
     &      (DMU**M5OVR3)*DMUP*DMUT-OVR3*U_NUC*(DMU**M2OVR3)*DMUPT      
                                                                        
          TMP3PX = OVR3*(DMU**M2OVR3)*DMUX-OVR3*M2OVR3*U_NUC*           &
     &      (DMU**M5OVR3)*DMUP*DMUX-OVR3*U_NUC*(DMU**M2OVR3)*DMUPX      
!                                                                       
!                                                                       
!                 Combination D function                                
!                                                                       
          SCRDU = U_NUC*COMPU*(TMP2+TMP3)/TMP1 
          SCRDUT = U_NUC*COMPU*(TMP2T+TMP3T)/TMP1 
          SCRDUX = U_NUC*COMPU*(TMP2X+TMP3X)/TMP1 
          SCRDXX = U_NUC*COMPU*(TMP2XX+TMP3XX)/TMP1 
          SCRDTT = U_NUC*COMPU*(TMP2TT+TMP3TT)/TMP1 
          SCRDXT = U_NUC*COMPU*(TMP2XT+TMP3XT)/TMP1 
!                                                                       
          SCRD = SCRDU/U_NUC 
          SCRDT = SCRDUT/U_NUC 
          SCRDX = SCRDUX/U_NUC 
!                                                                       
          SCRD2 = SCRDU/COMPU 
          SCRD2T = SCRDUT/COMPU 
          SCRD2X = SCRDUX/COMPU 
!                                                                       
          SCRDUP = SCRD-SCRD2+U_NUC*COMPU*                              &
     &    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)                
!                                                                       
          SCRDPT = SCRDT-SCRD2T+U_NUC*COMPU*                            &
     &    ((TMP2PT+TMP3PT)/TMP1-(TMP2T+TMP3T)*TMP1P/TMP1**2)            
!                                                                       
          SCRDPX = SCRDX-SCRD2X+U_NUC*COMPU*                            &
     &    ((TMP2PX+TMP3PX)/TMP1-(TMP2X+TMP3X)*TMP1P/TMP1**2)            
!                                                                       
          SCRDPP = (SCRDUP-SCRD)/U_NUC-(SCRD2+SCRDUP)/COMPU+            &
     &    (1.0-2.0*U_NUC)*                                              &
     &    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)+U_NUC*COMPU*   &
     &    ((TMP2PP+TMP3PP)/TMP1-2.0*(TMP2P+TMP3P)*TMP1P/TMP1**2-        &
     &    (TMP2+TMP3)*TMP1PP/TMP1**2+                                   &
     &    2.0*(TMP2+TMP3)*(TMP1P**2)/TMP1**3)                           
!                                                                       
!                                                                       
!                                                                       
!           bubble D function                                           
!bub          scrdu = (1.0-u_nuc)*dmu**ovr3                             
!bub          scrd = scrdu/u_nuc                                        
!bub          scrd2 = dmu**ovr3                                         
!bub          scrdup = -1.0*dmu**ovr3-                                  
!bub     1    ovr3*(1.0-u_nuc)*dmup*dmu**m2ovr3                         
!bub          scrdpp = ovr23*dmup*dmu**m2ovr3-ovr29*(1.0-u_nuc)*        
!bub     1    dmu**m5ovr3*dmup**2+ovr3*(1.0-u_nuc)*dmu**m2ovr3*dmupp    
!                                                                       
!                                                                       
!           nuclei D function                                           
!nuc          scrdu = u_nuc*du**ovr3                                    
!nuc          scrd = du**ovr3                                           
!nuc          scrd2 = scrdu/(1.0-u_nuc)                                 
!nuc          scrdup = du**ovr3+ovr3*u_nuc*dup*du**m2ovr3               
!nuc          scrdpp = ovr23*dup*du**m2ovr3-ovr29*u_nuc*                
!nuc     1    (du**m5ovr3)*(dup**2)+ovr3*u_nuc*(du**m2ovr3)*dupp        
!                                                                       
!                                                                       
!                                                                       
          ZETA_0 = CSSCAL*6.035204*(SIG_S*(16.0+Q))**OVR23 
!                                                                       
!                        Surface energy coefficent                      
          ZETA = ZETA_0*(H*SIGMA*X*NSUBI)**OVR23 
!                                                                       
!                        Derivative of Zeta w.r.t. X                    
          DZDT = OVR23*ZETA*HPRIM/H 
!                                                                       
!                        Derivative of Zeta w.r.t. X                    
          DZDX = OVR23*ZETA*(DHDX/H+SIGSGP+1.0/X) 
!                                                                       
!                        Derivative of Zeta w.r.t. NSUBI                
          DZDNI = OVR23*ZETA/NSUBI 
!                                                                       
!                                                                       
!                                                                       
!                        Nuclear radius                                 
          RSUBN = 9.0*H*SIGMA*SIG_0*U_NUC*(1.0-U_NUC)/                  &
     &    (2.0*ZETA*SCRDU)                                              
!                                                                       
!                        Nuclear volume                                 
          VSUBN = 4.0*PI*(RSUBN**3)/3.0 
!                                                                       
!                        Atomic number                                  
          A = NSUBI*VSUBN 
!                                                                       
!                        Now calc surface, Coulomb free energies        
!                                                                       
          FSUBSC = ZETA*SCRDU/BRYDNS 
          FSUBS = OVR23*ZETA*SCRDU/BRYDNS 
          FSUBC = OVR3*ZETA*SCRDU/BRYDNS 
!                                                                       
!                                                                       
!                                                                       
!                   Translational chemical potential                    
          MUSUBT = TRSCAL*                                              &
     &        T*DLOG((1.0-U_NUC)*(U_NUC*NSUBI)/(NQ*AZERO**2.5))         
!                                                                       
!                   Derivative of trans. chem. potential w.r.t. T       
          DMUTDT = TRSCAL*(MUSUBT/T-1.5) 
!                                                                       
!                   Translational free energy per baryon                
          FTRANS = TRSCAL*H*(MUSUBT-T)/AZERO 
!                                                                       
!                            if T is above the critical temperature     
        ELSE 
          A = 0.0 
          RSUBN = 0.0 
          VSUBN = 0.0 
          FSUBS = 0.0 
          FSUBC = 0.0 
          FTRANS = 0.0 
        ENDIF 
!                            Calc ratio of NSUBI to NSUBS               
        NRATIO = NSUBI/NSUBS 
!                                                                       
!                                                                       
        VNI = 2.0*AA*NSUBI+4.0*BB*X*NSUBI+CC*(1.0+DD)*NSUBI**DD 
!                                                                       
        VPI = 2.0*AA*NSUBI+4.0*BB*(1.0-X)*NSUBI+                        &
     &    CC*(1.0+DD)*NSUBI**DD+DELTAM                                  
!                                                                       
!---------------------------------------------------                    
!                                                                       
        ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ 
!                                                                       
        ZPI = 2.0*(PI**2)*NSUBI*X/MQ 
!                                                                       
        ETA_NI = FINV12(ZNI) 
!                                                                       
        ETA_PI = FINV12(ZPI) 
!                                                                       
        MUN_I = T*ETA_NI+VNI 
!                                                                       
        MUP_I = T*ETA_PI+VPI 
!                                                                       
        F32_NI = F_3_2(ETA_NI) 
!                                                                       
        F32_PI = F_3_2(ETA_PI) 
!                                                                       
        PSUBI = LQ*(F32_NI+F32_PI)+                                     &
     &    (NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+DD*CC*NSUBI**(1.0+DD)        
!                                                                       
!                                                                       
        BN = OVR23*ZETA*SCRD*(SIGSGP+HX+1.5*SCRDUX/SCRDU)*X/NSUBI-      &
     &  TRSCAL*(1.0-U_NUC)*(MUSUBT*(H-X*DHDX)/AZERO+X*DHDX*T/AZERO)     
!                                                                       
        BP = -OVR23*ZETA*SCRD*                                          &
     & ((SIGSGP+HX+1.5*SCRDUX/SCRDU)*COMPX+1.0/X)/NSUBI-                &
     & TRSCAL*(1.0-U_NUC)*                                              &
     & (MUSUBT*(H+DHDX*COMPX)/AZERO-DHDX*T*COMPX/AZERO)                 
!                                                                       
        BSUBP = ZETA*SCRDUP-OVR23*ZETA*SCRD-                            &
     &        TRSCAL*U_NUC*NSUBI*H*MUSUBT/AZERO                         
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!c        GPI = 2.0*FHALFI(ETA_PI)                                      
!c        GPI = 2.0*FHALFI(2.0*(pi**2)*x*nsubi/mq)                      
!c        GNI = 2.0*FHALFI(ETA_NI)                                      
!c        GNI = 2.0*FHALFI(2.0*(pi**2)*(1.0-x)*nsubi/mq)                
!                                                                       
!c        GPO = 2.0*FHALFO(ETA_PO)                                      
!c        GNO = 2.0*FHALFO(ETA_NO)                                      
!                                                                       
!                                                                       
        GPO = 2.0*FHALF(ETA_PO) 
        GNO = 2.0*FHALF(ETA_NO) 
        GPI = 2.0*FHALF(ETA_PI) 
        GNI = 2.0*FHALF(ETA_NI) 
!                                                                       
!                  Derivatives of inside potentials                     
!                                                                       
        DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
        DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
        DVNIDP = DVPIDN 
        DVNIDN = DVPIDP 
!                                                                       
!                  Derivatives of outside potentials                    
!                                                                       
        DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
        DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0))) 
        DVNODP = DVPODN 
        DVNODN = DVPODP 
!                                                                       
!                  Derivatives of inside K.E. densities                 
!                                                                       
        MSSCON = 3.0*MASSN/((HBAR*C)**2) 
        DTPIDP = MSSCON*T*GPI 
        DTPIDN = 0.0 
        DTNIDP = 0.0 
        DTNIDN = MSSCON*T*GNI 
!                                                                       
!                  Derivatives of outside K.E. densities                
!                                                                       
        DTPODP = MSSCON*T*GPO 
        DTPODN = 0.0 
        DTNODP = 0.0 
        DTNODN = MSSCON*T*GNO 
!                                                                       
!                                                                       
!                  Derivatives of inside chem. potentials               
!                                                                       
        DMPIDP = T*GPI/(X*NSUBI)+DVPIDP 
        DMPIDN = DVPIDN 
        DMNIDP = DVNIDP 
        DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN 
!                                                                       
!                  Derivatives of outside chem. potentials              
!                                                                       
        DMPODP = T+DVPODP*NPOUT/GPO 
        DMPODN = DVPODN*NNOUT/GNO 
        DMNODP = DVNODP*NPOUT/GPO 
        DMNODN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of inside pressure                       
!                                                                       
        DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP 
        DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN 
!                                                                       
!                  Derivatives of outside pressure                      
!                                                                       
        DPODP = NPOUT*DMPODP+NNOUT*DMNODP 
        DPODN = NPOUT*DMPODN+NNOUT*DMNODN 
!                                                                       
!                  Derivatives of alpha pressure                        
!                                                                       
        DPADP = ALFDNS*                                                 &
     &  ( (2.0-NPOUT*V_ALFA)*DMPODP+(2.0-NNOUT*V_ALFA)*DMNODP )         
        DPADN = ALFDNS*                                                 &
     &  ( (2.0-NPOUT*V_ALFA)*DMPODN+(2.0-NNOUT*V_ALFA)*DMNODN )         
!                                                                       
!                                                                       
        N1 = NSUBI-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
        N2 = NSUBI*X-EXALFA*NPOUT-2.0*ALFDNS 
!                                                                       
!                  Derivatives of U                                     
!                                                                       
        DUDPO = -EXCLU*(EXALFA*NPOUT/GPO+                               &
     &           (4.0-NOUT*V_ALFA)*DPADP/T)/N1                          
        DUDNO = -EXCLU*(EXALFA*NNOUT/GNO+                               &
     &           (4.0-NOUT*V_ALFA)*DPADN/T)/N1                          
        DUDNI = -U_NUC/N1 
!                                                                       
!                  Derivatives of X                                     
!                                                                       
        DXDPO = -(N2*DUDPO+EXCLU*(EXALFA*NPOUT/GPO+                     &
     &           (2.0-NPOUT*V_ALFA)*DPADP/T))/(U_NUC*NSUBI)             
        DXDNO = -(N2*DUDNO+EXCLU*(2.0-NPOUT*V_ALFA)*DPADN/T)/           &
     &           (U_NUC*NSUBI)                                          
        DXDNI = (N2-X*N1)/(NSUBI*N1) 
!                                                                       
!                  Derivatives of B's w.r.t. NSUBI                      
!                                                                       
        DB1DNI = TRSCAL*( -U_NUC*H*(MUSUBT+T)/AZERO )+                  &
     &      OVR23*ZETA*(SCRDUP-OVR23*SCRD)/NSUBI                        
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
!                                                                       
        DB2DNI = -2.0*ZETA*SCRD*TMP4/(9.0*NSUBI**2)-                    &
     &  TRSCAL*( (COMPU*T/(AZERO*NSUBI))*(H+COMPX*DHDX) )               
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU) 
        DB3DNI = -2.0*ZETA*SCRD*X*TMP4/(9.0*NSUBI**2)-                  &
     &          TRSCAL*( ((COMPU*T)/(AZERO*NSUBI))*(H-X*DHDX) )         
!                                                                       
!                                                                       
!                                                                       
!                  Derivatives of B's w.r.t. X                          
!                                                                       
        DB1DX = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*(SIGSGP+DHDX/H+1.0/X)+   &
     &  OVR23*ZETA*(SCRDPX-OVR23*SCRDX)-                                &
     &  TRSCAL*( U_NUC*NSUBI*DHDX*MUSUBT/AZERO )                        
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
!                                                                       
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+(X**(-2))+(X-1.0)*        &
     &  (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+1.5*SCRDXX/SCRDU-       &
     &  1.5*(SCRDUX/SCRDU)**2)                                          
!                                                                       
!                                                                       
        DB2DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+      &
     &      OVR23*ZETA*SCRD*TMP5/NSUBI-TRSCAL*(                         &
     &      COMPU*(DHDX*MUSUBT+(DHDXX*(1.0-X)-DHDX)*(MUSUBT-T))/AZERO)  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X 
!                                                                       
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+X*                        &
     &         (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+                 &
     &       1.5*SCRDXX/SCRDU-1.5*(SCRDUX/SCRDU)**2)                    
!                                                                       
        DB3DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+      &
     &      OVR23*ZETA*SCRD*TMP5/NSUBI-                                 &
     &      TRSCAL*( COMPU*(DHDX*T-X*DHDXX*(MUSUBT-T))/AZERO )          
!                                                                       
!                                                                       
!                                                                       
!                  Derivatives of B's w.r.t. U_NUC                      
!                                                                       
        DB1DU = ZETA*(SCRDPP-OVR23*SCRDUP/U_NUC+OVR23*SCRD/U_NUC)-      &
     &  TRSCAL*( NSUBI*H*(MUSUBT+T*(1.0-2.0*U_NUC)/(1.0-U_NUC))/AZERO ) 
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
        TMP5 = (X-1.0)*1.5*(SCRDPX/SCRDU-SCRDUX*SCRDUP/SCRDU**2) 
        DB2DU = (OVR23*ZETA*SCRD/NSUBI)*TMP4*(SCRDUP/SCRDU-1.0/U_NUC)+  &
     &    OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+                          &
     &    TRSCAL*( (H*MUSUBT+DHDX*COMPX*(MUSUBT-T))/AZERO-              &
     &    (T*(1.0-2.0*U_NUC)/U_NUC)*(H+DHDX*COMPX)/AZERO )              
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X 
        TMP5 = X*1.5*(SCRDPX/SCRDU-SCRDUP*SCRDUX/SCRDU**2) 
        DB3DU = OVR23*ZETA*SCRD*TMP4*(U_NUC*SCRDUP/SCRDU-1.0)/          &
     & (U_NUC*NSUBI)+OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+               &
     &  TRSCAL*( (H*MUSUBT-X*DHDX*(MUSUBT-T))/AZERO-                    &
     & T*(1.0-2.0*U_NUC)*(H-X*DHDX)/(AZERO*U_NUC) )                     
!                                                                       
!                                                                       
!                      A1 derivatives                                   
!                                                                       
        DA1ID1 = X*DPIDP+(1.0-X)*DPIDN+NSUBI*(DPIDP-DPIDN)*DXDNI 
        DA1ID2 = NSUBI*(DPIDP-DPIDN)*DXDPO 
        DA1ID3 = NSUBI*(DPIDP-DPIDN)*DXDNO 
!                                                                       
        DA1OD1 = 0.0 
        DA1OD2 = DPODP+DPADP 
        DA1OD3 = DPODN+DPADN 
!                                                                       
        DB1D1 = DB1DNI+DB1DX*DXDNI+DB1DU*DUDNI 
        DB1D2 = DB1DX*DXDPO+DB1DU*DUDPO 
        DB1D3 = DB1DX*DXDNO+DB1DU*DUDNO 
!                                                                       
        DA1D1 = DA1ID1-DB1D1-DA1OD1 
        DA1D2 = DA1ID2-DB1D2-DA1OD2 
        DA1D3 = DA1ID3-DB1D3-DA1OD3 
!                                                                       
!                      A3 derivatives                                   
!                                                                       
        DA3ID1 = X*DMNIDP+(1.0-X)*DMNIDN+NSUBI*(DMNIDP-DMNIDN)*DXDNI 
        DA3ID2 = NSUBI*(DMNIDP-DMNIDN)*DXDPO 
        DA3ID3 = NSUBI*(DMNIDP-DMNIDN)*DXDNO 
!                                                                       
        DA3OD1 = 0.0 
        DA3OD2 = DMNODP 
        DA3OD3 = DMNODN 
!                                                                       
        DB3D1 = DB3DNI+DB3DX*DXDNI+DB3DU*DUDNI 
        DB3D2 = DB3DX*DXDPO+DB3DU*DUDPO 
        DB3D3 = DB3DX*DXDNO+DB3DU*DUDNO 
!                                                                       
        DA3D1 = DA3ID1-DB3D1-DA3OD1 
        DA3D2 = DA3ID2-DB3D2-DA3OD2 
        DA3D3 = DA3ID3-DB3D3-DA3OD3 
!                                                                       
!                      A2 derivatives                                   
!                                                                       
        DA2ID1 = X*DMPIDP+(1.0-X)*DMPIDN+NSUBI*(DMPIDP-DMPIDN)*DXDNI 
        DA2ID2 = NSUBI*(DMPIDP-DMPIDN)*DXDPO 
        DA2ID3 = NSUBI*(DMPIDP-DMPIDN)*DXDNO 
!                                                                       
        DA2OD1 = 0.0 
        DA2OD2 = DMPODP 
        DA2OD3 = DMPODN 
!                                                                       
        DB2D1 = DB2DNI+DB2DX*DXDNI+DB2DU*DUDNI 
        DB2D2 = DB2DX*DXDPO+DB2DU*DUDPO 
        DB2D3 = DB2DX*DXDNO+DB2DU*DUDNO 
!                                                                       
        DA2D1 = DA2ID1-DB2D1-DA2OD1 
        DA2D2 = DA2ID2-DB2D2-DA2OD2 
        DA2D3 = DA2ID3-DB2D3-DA2OD3 
!                                                                       
!                                                                       
!                      Eta derivatives                                  
!                                                                       
        DNDETN = NNOUT/GNO 
        DPDETP = NPOUT/GPO 
!                                                                       
        DA1DN = DA1D1 
        DA1ETP = DA1D2 
        DA1ETN = DA1D3 
!                                                                       
        DA2DN = DA2D1 
        DA2ETP = DA2D2 
        DA2ETN = DA2D3 
!                                                                       
        DA3DN = DA3D1 
        DA3ETP = DA3D2 
        DA3ETN = DA3D3 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        A1 = PSUBI-BSUBP-BPROUT-BPRALF 
        A2 = MUP_I-BP-MUP_O 
        A3 = MUN_I-BN-MUN_O 
!                                                                       
!                                                                       
!                          Unset the "new" flag                         
        NEWFLG = 0 
!                                                                       
        DETERM = DA1DN*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)-                   &
     &           DA1ETP*(DA2DN*DA3ETN-DA2ETN*DA3DN)+                    &
     &           DA1ETN*(DA2DN*DA3ETP-DA2ETP*DA3DN)                     
!                                                                       
        DNSUBI = -1.0*(A1*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)+                &
     &           A2*(DA3ETP*DA1ETN-DA1ETP*DA3ETN)+                      &
     &           A3*(DA1ETP*DA2ETN-DA1ETN*DA2ETP))/DETERM               
!                                                                       
!                                                                       
        DETAP = -1.0*(A1*(DA2ETN*DA3DN-DA2DN*DA3ETN)+                   &
     &          A2*(DA1DN*DA3ETN-DA1ETN*DA3DN)+                         &
     &          A3*(DA1ETN*DA2DN-DA1DN*DA2ETN))/DETERM                  
!                                                                       
!                                                                       
        DETAN = -1.0*(A1*(DA2DN*DA3ETP-DA2ETP*DA3DN)+                   &
     &          A2*(DA1ETP*DA3DN-DA1DN*DA3ETP)+                         &
     &          A3*(DA1DN*DA2ETP-DA1ETP*DA2DN))/DETERM                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                        Check the step size in NSUBI                   
        IF(ABS(DNSUBI/NSUBI).GT.0.04) THEN 
          DNSUBI = 0.04*DNSUBI*NSUBI/ABS(DNSUBI) 
        ENDIF 
   26   CONTINUE 
        NSUBIN = NSUBI+DNSUBI 
        IF((NSUBIN.LT.DMAX1(4.5D-2,BRYDNS)).OR.(NSUBIN.GT.0.25)) THEN 
          DNSUBI = 0.5*DNSUBI 
          GOTO 26 
        ENDIF 
!                                                                       
!                        Check the step size in ETA_PO                  
        IF(ABS(DETAP).GT.4.0) THEN 
          DETAP = 4.0*DETAP/ABS(DETAP) 
        ENDIF 
   27   CONTINUE 
        NETAP = ETA_PO+DETAP 
        IF((NETAP.LT.-5000.0).OR.(NETAP.GT.ETAMAX)) THEN 
          DETAP = 0.5*DETAP 
          GOTO 27 
        ENDIF 
!                                                                       
!                        Check the step size in ETA_NO                  
        IF(ABS(DETAN).GT.4.0) THEN 
          DETAN = 4.0*DETAN/ABS(DETAN) 
        ENDIF 
   28   CONTINUE 
        NETAN = ETA_NO+DETAN 
        IF((NETAN.LT.-5000.0).OR.(NETAN.GT.ETAMAX)) THEN 
          DETAN = 0.5*DETAN 
          GOTO 28 
        ENDIF 
!                                                                       
!                                                                       
!                        Update the variables                           
!cc        if(i.lt.30) write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc      
 1205   format(i3,1p9e21.14) 
!                                                                       
        NSUBI = NSUBI+DNSUBI 
        ETA_PO = ETA_PO+DETAP 
        ETA_NO = ETA_NO+DETAN 
!                                                                       
!                                                                       
!                                                                       
!                        If the required tolarences have been met       
!                        break out of the loop                          
        IF((ABS(DNSUBI).LT.NSIACC).AND.(ABS(DETAP).LT.PRTACC)           &
     &    .AND.(ABS(DETAN).LT.NUTACC) ) THEN                            
          GOTO 40 
        ELSE 
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '2',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
          GOTO 30 
        ENDIF 
!                                                                       
!                                                                       
   29   CONTINUE 
        IF(NEWFLG.NE.1) THEN 
          cflag = cflag+1 
          DNSUBI = 0.5*DNSUBI 
          NSUBI = NSUBI-DNSUBI 
          DETAP = 0.5*DETAP 
          ETA_PO = ETA_PO-DETAP 
          DETAN = 0.5*DETAN 
          ETA_NO = ETA_NO-DETAN 
          IF(DBFLAG.EQ.1) THEN 
            WRITE(*,2000) '3',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
          ENDIF 
          GOTO 30 
        ELSE 
          NSUBI = NSUBS 
!c          ETA_PO = ETA_PO-0.5/T                                       
!c          ETA_NO = ETA_NO-0.5/T                                       
          ETA_PO = ETA_PO-2.0/T 
          ETA_NO = ETA_NO-2.0/T 
        ENDIF 
!                                                                       
!                                                                       
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '4',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
 2000   FORMAT(t2,a,1x,i3,1x,f8.5,3(1X,G13.5)) 
!                                                                       
!                                                                       
   30 END DO 
!                                                                       
!            If scheme 1 has failed try scheme 2                        
      if(schflg.eq.0) then 
        schflg = 1 
        goto 5 
      endif 
!                                                                       
!                                                                       
      SSFLAG = 0 
      GOTO 999 
!                                                                       
!                    Branch label to break out of DO 30 iteration       
   40 CONTINUE 
!                                                                       
!                                                                       
!                    The following logic determines whether this was    
!                    the correct scheme to use, and if not then which   
!                    one should be used                                 
!                                                                       
      if(ftflag.ne.0) then 
        ssflag = 4 
        goto 999 
      endif 
!                                                                       
!                    If calculated critical temperature is less than T, 
!                    then switch to the scheme with no nuclei           
      IF(T.GE.TSUBC) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                    If fraction of nuclei present is zero and no switch
!                    has been made then switch to the no nuclei scheme  
      IF(U_NUC.LE.0.0) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ELSEIF(U_NUC.GT.1.0) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ELSE 
!                    Set flag to indicate success                       
        SSFLAG = 1 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                    If eqns aren't really zeroed then fail             
!                                                                       
!                                                                       
      IF( (ABS(A1).GT.1.0D-5).OR.(ABS(A2).GT.1.0D-5).OR.                &
     &    (ABS(A3).GT.1.0D-5) ) THEN                                    
        SSFLAG = 0 
!c        WRITE(*,*) ' NUCEOS: False convg; A = ',A1,A2,A3              
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF(NSUBI.LT.0.05) THEN 
        WRITE(*,*) 'NUCEOS:: <<WARNING>> NSUBI GETTING CLOSE TO LB' 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
      ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ 
!                                                                       
      ZPI = 2.0*(PI**2)*NSUBI*X/MQ 
!                                                                       
      ETA_NI = FINV12(ZNI) 
!                                                                       
      ETA_PI = FINV12(ZPI) 
!                                                                       
      MUN_I = T*ETA_NI+VNI 
!                                                                       
      MUP_I = T*ETA_PI+VPI 
!                                                                       
      F32_NI = F_3_2(ETA_NI) 
!                                                                       
      F32_PI = F_3_2(ETA_PI) 
!                                                                       
      EXCLU = 1.0-U_NUC 
      EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                                                                       
!                                                                       
!                    Calculate particle fractions                       
!                                                                       
      XALFA = 4.0*EXCLU*ALFDNS/BRYDNS 
      XNUT = NNOUT*EXCLU*EXALFA/BRYDNS 
      XPROT = NPOUT*EXCLU*EXALFA/BRYDNS 
      XH = 1.0-XPROT-XNUT-XALFA 
      XHCHK = U_NUC*NSUBI/BRYDNS 
!                                                                       
      IF((XH.LT.HEAVCT).OR.(XHCHK.LT.HEAVCT)) THEN 
!                    Set flag to indicate switch is being made          
        SSFLAG = 0 
!c        write(*,*) ' xh,xhchk = ',xh,xhchk                            
        GOTO 999 
      ENDIF 
!                                                                       
      IF((XALFA.LT.0.0).OR.(XH.LT.0.0).OR.                              &
     &   (XNUT.LT.0.0).OR.(XPROT.LT.0.0)) THEN                          
        SSFLAG = 0 
        write(*,*) ' Xs hnpa = ',xh,xnut,xprot,xalfa 
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Baryons                                            
!                                                                       
!                                                                       
      MUPROT = MUP_O 
      MUN = MUN_O 
      MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
      IF(ABS((XH-XHCHK)/XHCHK).GT.1.0D-4) THEN 
        SSFLAG = 0 
        GOTO 999 
!CC        WRITE(*,*) ' INCONSISTENCEY IN XH AT',T,BRYDNS,YE,XH,XHCHK   
      ENDIF 
!                                                                       
      NUCDNS = BRYDNS*XH 
!                                                                       
      TAU_PO = KQ*F32_PO 
      TAU_PI = KQ*F32_PI 
!                                                                       
      TAU_NO = KQ*F32_NO 
      TAU_NI = KQ*F32_NI 
!                                                                       
      IF(NOUT.GT.0.0) XOUT = NPOUT/NOUT 
!                                                                       
!                                                                       
!                    Calculate internal energy of outside nucleons,     
!                    alpha particles, and nuclei (per baryon)           
!                                                                       
      BUOUT = (EXCLU*EXALFA/BRYDNS)*( UQ*(TAU_PO+TAU_NO)+EIFLAG*        &
     &    ( (NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+                            &
     &    CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )                             
!                                                                       
      BUNUC = XH*( ( UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*                     &
     & (AA+4.0*BB*X*(1.0-X))+CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM )/       &
     & NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+              &
     & TRSCAL*                                                          &
     & (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)   
!                                                                       
!                                                                       
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA) 
!                                                                       
      BU = BUOUT+BUALFA+BUNUC 
!                                                                       
!                                                                       
      BSOUT = (EXCLU*EXALFA/BRYDNS)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)- &
     & NNOUT*ETA_NO-NPOUT*ETA_PO )                                      
!                                                                       
!                                                                       
!                    Calculate entropy of alpha particles (per baryon)  
      BSALFA = -0.25*XALFA*(MUALFA/T-2.5) 
!                                                                       
!                                                                       
      BSNUC = XH*( (5.0*UQ/(3.0*T))*(TAU_NI+TAU_PI)-                    &
     & NSUBI*(1.0-X)*ETA_NI-NSUBI*X*ETA_PI )/NSUBI-                     &
     & FSUBSC*(SCRDUT/SCRDU+OVR23*HPRIM/H)-                             &
     & XH*TRSCAL*(1.0-U_NUC)*                                           &
     & ((FTRANS*HPRIM/H)+H*(MUSUBT/T-2.5)/AZERO)                        
!                                                                       
!                    Calculate total baryon entropy (per baryon)        
      BS = BSOUT+BSNUC+BSALFA 
!                                                                       
!                    Calculate free energy of outside nucleons (per bary
      BFOUT = BUOUT-T*BSOUT 
!                                                                       
!                    Calculate free energy of alpha particles (per baryo
      BFALFA = BUALFA-T*BSALFA 
!                                                                       
!                    Calculate free energy of nuclei (per baryon)       
      BFNUC = BUNUC-T*BSNUC 
!                                                                       
!                    Calculate total baryon free energy (per baryon)    
      BFTOT = BFOUT+BFNUC+BFALFA 
!                                                                       
!                    Calculate pressure due to nuclei                   
      BPRNUC = -ZETA*(SCRDU-U_NUC*SCRDUP)+                              &
     & TRSCAL*U_NUC*NSUBI*H*((1.0-U_NUC)*T-U_NUC*MUSUBT)/AZERO          
!                                                                       
!                                                                       
!                    Calculate total baryon pressure                    
      BPRESS = BPROUT+BPRALF+BPRNUC 
!                                                                       
!                                                                       
!                    Leptons & Photons                                  
!                                                                       
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
!                                                                       
!                    Total pressure and eng/ent per baryon              
!                                                                       
      FBARY = BFTOT+FSUBE 
      PBARY = BPRESS+EPRESS 
      MUBARY = YE*MUPROT+(1.0-YE)*MUN 
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN 
!                                                                       
      FTOT = BFTOT+FSUBE+PF 
      UTOT = BU+EU+PU 
      STOT = BS+ES+PS 
      PTOT = BPRESS+EPRESS+PPRESS 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
!                Derivatives of thermodynamic variables                 
!-----------------------------------------------------------------------
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of exterior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and ETA's)    !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                  Derivatives of exterior potentials                   
!                  w.r.t. particle densities                            
      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVNODP = DVPODN 
      DVNODN = DVPODP 
!                                                                       
!                                                                       
!                  Derviatives of exterior chem. pot. w.r.t. ETA's      
!                  (at fixed T)                                         
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of pressure potential w.r.t.             
!                  particle densities                                   
      DV_DPO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
      DV_DNO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
!                                                                       
!                  Derivatives of pressure potential w.r.t. ETA's       
!                  (at fixed T)                                         
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
!                  Derivatives of outside pressure w.r.t. ETA's         
!                  (at fixed T)                                         
      DPODEP = NPOUT*T+DV_DEP 
      DPODEN = NNOUT*T+DV_DEN 
!                                                                       
!                  Derivatives of alpha density w.r.t. ETA's            
!                  (at fixed T)                                         
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T 
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T 
!                                                                       
!                  Derivatives of alpha pressure w.r.t. ETA's           
!                  (at fixed T)                                         
      DPADEP = T*DNADEP 
      DPADEN = T*DNADEN 
!                                                                       
!                  Derivatives of particle densities w.r.t. T           
!                  (at fixed ETA's)                                     
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
!                  Derivatives of exterior chem. pot. w.r.t. T          
!                  (at fixed ETA's)                                     
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
!                  (at fixed ETA's)                                     
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
!                  (at fixed ETA's)                                     
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                  Derivative of alpha chem. pot. w.r.t. T              
!                  (at fixed ETA's)                                     
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
!                  (at fixed ETA's)                                     
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                  Derivative of alpha particle pressure w.r.t. T       
!                  (at fixed ETA's)                                     
      DPADT = ALFDNS+T*DNADT 
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of interior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and density)  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                   Derivatives of kinetic energy densities w.r.t. T    
!                   (holding the number densities (X & NSUBI) fixed)    
      DTPIDT =2.5*TAU_PI/T-2.25*X*NSUBI*GPI/UQ 
      DTNIDT =2.5*TAU_NI/T-2.25*(1.0-X)*NSUBI*GNI/UQ 
!                                                                       
!                   Derivatives of pressures w.r.t. T                   
!                   (holding the number densities (X & NSUBI) fixed)    
      DPIDT = OVR23*UQ*(DTPIDT+DTNIDT) 
!                                                                       
!                   Derivatives of interior chem. pot. w.r.t. T         
!                   (holding the number densities (X & NSUBI) fixed)    
      DMPIDT = ETA_PI-1.5*GPI 
      DMNIDT = ETA_NI-1.5*GNI 
!                                                                       
!                                                                       
!                  Derivatives of inside potentials w.r.t.              
!                  interior proton and neutron densities                
!                  (at fixed T)                                         
      DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
      DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
      DVNIDP = DVPIDN 
      DVNIDN = DVPIDP 
!                                                                       
!                   Derivatives of interior chemical potentials         
!                   w.r.t. interior neutron and proton densities        
!                  (at fixed T)                                         
      DMPIDP = T*GPI/(X*NSUBI)+DVPIDP 
      DMPIDN = DVPIDN 
      DMNIDP = DVNIDP 
      DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN 
!                                                                       
!                   Derivatives of interior pressure                    
!                   w.r.t. interior neutron and proton densities        
!                  (at fixed T)                                         
      DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP 
      DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of "B" terms    !                  
!                 !      from the chemical and       !                  
!                 !      pressure equilibrium        !                  
!                 !      equations                   !                  
!                 !                                  !                  
!                 !      (w.r.t. Temperature )       !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!             Derivative of term from pressure equilibrium eqn.         
!                                                                       
      DB1DT = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*HPRIM/H+                   &
     &    ZETA*(SCRDPT-OVR23*SCRDT)-                                    &
     &    TRSCAL*U_NUC*NSUBI*(HPRIM*MUSUBT+H*DMUTDT)/AZERO              
!                                                                       
!                                                                       
!             Derivative of term from proton equilibrium eqn.           
!                                                                       
      TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+                                  &
     & 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2                      
!                                                                       
      DB2DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*TMP4+                   &
     &    OVR23*ZETA*SCRDT*TMP4/NSUBI+                                  &
     &    OVR23*(ZETA*SCRD/NSUBI)*(X-1.0)*TMP5-                         &
     &    TRSCAL*EXCLU*(DMUTDT*(H+DHDX*(1.0-X))+MUSUBT*                 &
     &    (HPRIM+DHDTDX*(1.0-X))-DHDX*(1.0-X)-T*DHDX*(1.0-X))/AZERO     
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!             Derivative of term from neutron equilibrium eqn.          
!                                                                       
      TMP4 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU 
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+                                  &
     & 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2                      
      DB3DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*X*TMP4+                 &
     &        OVR23*(ZETA*SCRDT/NSUBI)*X*TMP4+                          &
     &        OVR23*(ZETA*SCRD/NSUBI)*X*TMP5-                           &
     &        TRSCAL*EXCLU*(HPRIM*MUSUBT+H*DMUTDT-X*DHDTDX*(MUSUBT-T)-  &
     &        X*DHDX*(DMUTDT-1.0))/AZERO                                
!                                                                       
!                                                                       
!                                                                       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of constraint   !                  
!                 !      and equilibrium equations   !                  
!                 !      with respect to the five    !                  
!                 !      compositional variables     !                  
!                 !      (U,x,n_i,eta_po,eta_no)     !                  
!                 !      and the three independent   !                  
!                 !      variables                   !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 1 (Baryon conservation)                       
!                                                                       
      DFDOM(1,1) = NOUT*EXALFA+4.0*ALFDNS-NSUBI 
!                                                                       
      DFDOM(1,2) = 0.0 
!                                                                       
      DFDOM(1,3) = -U_NUC 
!                                                                       
      DFDOM(1,4) = -EXCLU*EXALFA*NPOUT/GPO+                             &
     &             V_ALFA*DNADEP*EXCLU*NOUT-4.0*EXCLU*DNADEP            
!                                                                       
      DFDOM(1,5) = -EXCLU*EXALFA*NNOUT/GNO+                             &
     &             V_ALFA*DNADEN*EXCLU*NOUT-4.0*EXCLU*DNADEN            
!                                                                       
!                                                                       
!                                                                       
      DFDL_1(1) = -1.0 
!                                                                       
      DFDL_2(1) = EXCLU*EXALFA*(DNPODT+DNNODT)-EXCLU*V_ALFA*NOUT*DNADT+ &
     &     4.0*EXCLU*DNADT                                              
!                                                                       
      DFDL_3(1) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 2 (Charge conservation)                       
!                                                                       
      DFDOM(2,1) = EXALFA*NPOUT+2.0*ALFDNS-X*NSUBI 
!                                                                       
      DFDOM(2,2) = -U_NUC*NSUBI 
!                                                                       
      DFDOM(2,3) = -X*U_NUC 
!                                                                       
      DFDOM(2,4) = -EXCLU*EXALFA*NPOUT/GPO+                             &
     &     V_ALFA*EXCLU*NPOUT*DNADEP-2.0*EXCLU*DNADEP                   
!                                                                       
      DFDOM(2,5) = V_ALFA*EXCLU*NPOUT*DNADEN-2.0*EXCLU*DNADEN 
!                                                                       
!                                                                       
!                                                                       
      DFDL_1(2) = -1.0*YE 
!                                                                       
      DFDL_2(2) = EXCLU*EXALFA*DNPODT-V_ALFA*EXCLU*NPOUT*DNADT+         &
     &     2.0*EXCLU*DNADT                                              
!                                                                       
      DFDL_3(2) = -1.0*BRYDNS 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 3 (Proton chemical equilibrium)               
!                                                                       
      DFDOM(3,1) = -DB2DU 
!                                                                       
      DFDOM(3,2) = NSUBI*(DMPIDP-DMPIDN)-DB2DX 
!                                                                       
      DFDOM(3,3) = (1.0-X)*DMPIDN+X*DMPIDP-DB2DNI 
!                                                                       
      DFDOM(3,4) = -DMPDEP 
!                                                                       
      DFDOM(3,5) = -DMPDEN 
!                                                                       
      DFDL_1(3) = 0.0 
      DFDL_2(3) = -1.0*(DMPIDT-DMPODT-DB2DT) 
      DFDL_3(3) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 4 (Neutron chemical equilibrium)              
!                                                                       
      DFDOM(4,1) = -DB3DU 
!                                                                       
      DFDOM(4,2) = NSUBI*(DMNIDP-DMNIDN)-DB3DX 
!                                                                       
      DFDOM(4,3) = (1.0-X)*DMNIDN+X*DMNIDP-DB3DNI 
!                                                                       
      DFDOM(4,4) = -DMNDEP 
!                                                                       
      DFDOM(4,5) = -DMNDEN 
!                                                                       
      DFDL_1(4) = 0.0 
      DFDL_2(4) = -1.0*(DMNIDT-DMNODT-DB3DT) 
      DFDL_3(4) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 5 (Pressure equilibrium)                      
!                                                                       
      DFDOM(5,1) = -DB1DU 
!                                                                       
      DFDOM(5,2) = NSUBI*(DPIDP-DPIDN)-DB1DX 
!                                                                       
      DFDOM(5,3) = (1.0-X)*DPIDN+X*DPIDP-DB1DNI 
      ncomp = dfdom(5,3) 
!                                                                       
      DFDOM(5,4) = -DPODEP-DPADEP 
!                                                                       
      DFDOM(5,5) = -DPODEN-DPADEN 
!                                                                       
      DFDL_1(5) = 0.0 
      DFDL_2(5) = -1.0*(DPIDT-DPODT-DPADT-DB1DT) 
      DFDL_3(5) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!      write(*,*) ' '                                                   
!c      write(*,*) ' '                                                  
!c      write(*,7010) db1dx,db2dx,db3dx                                 
!c      write(*,7010) db1du,db2du,db3du                                 
!c      write(*,7010) db1dni,db2dni,db3dni                              
!c      write(*,7010) db1dt,db2dt,db3dt                                 
 7010 format(3(1x,g13.6)) 
!      write(*,7000) x,u_nuc,nsubi,eta_po,eta_no                        
!      write(*,*) 'eta_i ',eta_pi,eta_ni                                
!      write(*,*) ' as ',a1,a2,a3                                       
!      write(*,*) ' '                                                   
!      write(*,7000) (dfdom(1,i),i=1,5,1)                               
!      write(*,7000) (dfdom(2,i),i=1,5,1)                               
!      write(*,7000) (dfdom(3,i),i=1,5,1)                               
!      write(*,7000) (dfdom(4,i),i=1,5,1)                               
!      write(*,7000) (dfdom(5,i),i=1,5,1)                               
!      write(*,*) ' '                                                   
!c      write(*,*) ' dna: ',dnadpo,dnadno                               
!c      write(*,*) ' dt: ',dmpidt,dmpodt                                
!      write(*,7000) (dfdl_1(i),i=1,5,1)                                
 7000 format(5(1x,g13.6)) 
!                                                                       
!      pause                                                            
!                    Invert the DFDOM matrix                            
!                                                                       
      CALL MATINV(DFDOM,DFDOMI,5) 
!  IMSL subroutine call to invert the matrix                            
!CC      CALL DLINRG(5,DFDOM,5,DFDOMI,5)                                
!                                                                       
!c      call matmul(dfdom,dfdomi,a_tmp,5,5,5)                           
!                                                                       
!c      write(*,*) ' '                                                  
!c      write(*,7000) (a_tmp(1,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(2,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(3,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(4,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(5,i),i=1,5,1)                              
!                                                                       
!c      write(*,7000) (dfdomi(1,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(2,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(3,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(4,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(5,i),i=1,5,1)                             
!c      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '      
!c      write(*,7000) (dfdl_2(i),i=1,5,1)                               
!c      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '      
!                                                                       
!      DO 800 LLI=1,5,1                                                 
!        R_CHECK(LLI) = 0.0                                             
!        DO 801 KKI=1,5,1                                               
!          R_CHECK(LLI) = R_CHECK(LLI)+DFDOM(LLI,KKI)*RESULT(KKI)       
! 801    CONTINUE                                                       
!        r_check(lli) = r_check(lli)-dfdl_1(lli)                        
! 800  CONTINUE                                                         
!      write(*,*) ' >>>>>>>>>>>>>>>    R check <<<<<<<<<<<<<<<<< '      
!      write(*,7000) (r_check(i),i=1,5,1)                               
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                    Multiply the DFDL_1 vector by the DFDOMI matrix    
!                    to get the density derivatives                     
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_1,RESULT,5) 
!                                                                       
      DU_DN = RESULT(1) 
      DX_DN = RESULT(2) 
      DNI_DN = RESULT(3) 
      DEP_DN = RESULT(4) 
      DEN_DN = RESULT(5) 
!                                                                       
!                                                                       
!                    Multiply the DFDL_2 vector by the DFDOMI matrix    
!                    to get the Temperature derivatives                 
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_2,RESULT,5) 
!                                                                       
      DU_DT = RESULT(1) 
      DX_DT = RESULT(2) 
      DNI_DT = RESULT(3) 
      DEP_DT = RESULT(4) 
      DEN_DT = RESULT(5) 
!                                                                       
!                    Multiply the DFDL_3 vector by the DFDOMI matrix    
!                    to get the Ye derivatives                          
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_3,RESULT,5) 
!                                                                       
      DU_DY = RESULT(1) 
      DX_DY = RESULT(2) 
      DNI_DY = RESULT(3) 
      DEP_DY = RESULT(4) 
      DEN_DY = RESULT(5) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of finite size  !                  
!                 !      terms in the internal       !                  
!                 !      energy and entropy          !                  
!                 !      densities w.r.t. to U,X,n_i !                  
!                 !      and T.  These are used in   !                  
!                 !      calculating the derivatives !                  
!                 !      w.r.t. the independant vars !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                        Free energy Surface & Coulomb terms            
!                                  (Densities)                          
!                                                                       
      F_SC = ZETA*SCRDU 
!                                                                       
      DFSCDU = ZETA*SCRDUP 
!                                                                       
      DFSCDX = ZETA*SCRDUX+SCRDU*DZDX 
!                                                                       
      DFSCDN = SCRDU*DZDNI 
!                                                                       
      DFSCDT = ZETA*SCRDUT+SCRDU*DZDT 
!                                                                       
!                                                                       
!                        Free energy translational terms                
!                                  (Densities)                          
      FTR = U_NUC*EXCLU*NSUBI*FTRANS 
!                                                                       
      DFTRDT = FTR*(HPRIM/H+1.0/T)-                                     &
     &    1.5*TRSCAL*U_NUC*EXCLU*NSUBI*H/AZERO                          
!                                                                       
      DFTRDX = FTR*DHDX/H 
!                                                                       
      DFTRDU = FTR/U_NUC-FTR/EXCLU+                                     &
     &    TRSCAL*NSUBI*H*(1.0-2.0*U_NUC)/AZERO                          
!                                                                       
      DFTRDN = FTR/NSUBI+TRSCAL*U_NUC*EXCLU*H*T/AZERO 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!                                                                       
!                        Internal energy Surface & Coulomb terms        
!                                  (Densities)                          
!                                                                       
      TMP4 = 1.0-T*SCRDUT/SCRDU-OVR23*T*HPRIM/H 
!                                                                       
      E_SC = F_SC*TMP4 
!                                                                       
      DESCDU = DFSCDU*TMP4+                                             &
     &    F_SC*(T*SCRDUT*SCRDUP/SCRDU**2-T*SCRDPT/SCRDU)                
!                                                                       
      DESCDX = DFSCDX*TMP4+                                             &
     &    F_SC*(T*SCRDUT*SCRDUX/SCRDU**2-T*SCRDXT/SCRDU+                &
     &    OVR23*T*HPRIM*DHDX/H**2-OVR23*T*DHDTDX/H)                     
!                                                                       
      DESCDN = DFSCDN*TMP4 
!                                                                       
      DESCDT = DFSCDT*TMP4+F_SC*                                        &
     &   (T*(SCRDUT**2)/SCRDU**2-SCRDUT/SCRDU-T*SCRDTT/SCRDU+           &
     &    OVR23*T*(HPRIM**2)/H**2-OVR23*HPRIM/H-OVR23*T*HPPRIM/H)       
!                                                                       
!                        Internal energy translational terms            
!                                  (Densities)                          
!                                                                       
      TMP4 = 1.5*H*T/AZERO-T*HPRIM*(MUSUBT-T)/AZERO 
!                                                                       
      E_TR = TRSCAL*EXCLU*BRYDNS*XH*TMP4 
!                                                                       
      DETRDU = TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4-                      &
     &    NSUBI*(T**2)*HPRIM*(1.0-2.0*U_NUC)/AZERO)                     
!                                                                       
      DETRDX = TRSCAL*BRYDNS*XH*EXCLU*                                  &
     &    (1.5*T*DHDX/AZERO-T*(MUSUBT-T)*DHDTDX/AZERO)                  
!                                                                       
      DETRDN = TRSCAL*(U_NUC*EXCLU*TMP4-                                &
     &    BRYDNS*XH*EXCLU*(T**2)*HPRIM/(NSUBI*AZERO))                   
!                                                                       
      DETRDT = TRSCAL*BRYDNS*XH*EXCLU*                                  &
     &    (1.5*(H+T*HPRIM)/AZERO-(HPRIM+T*HPPRIM)*(MUSUBT-T)/AZERO-     &
     &    T*HPRIM*(MUSUBT/T-2.5)/AZERO )                                
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!                                                                       
!                        Entropy Surface & Coulomb terms                
!                                  (Densities)                          
!                                                                       
      S_SC = (E_SC-F_SC)/T 
!                                                                       
      DSSCDU = (DESCDU-DFSCDU)/T 
!                                                                       
      DSSCDX = (DESCDX-DFSCDX)/T 
!                                                                       
      DSSCDN = (DESCDN-DFSCDN)/T 
!                                                                       
      DSSCDT = (DESCDT-DFSCDT)/T-(E_SC-F_SC)/T**2 
!                                                                       
!                        Entropy translational terms                    
!                                  (Densities)                          
!                                                                       
      TMP4 = MUSUBT*(HPRIM+H/T)/AZERO-(T*HPRIM+2.5*H)/AZERO 
!                                                                       
      S_TR = -TRSCAL*BRYDNS*XH*EXCLU*TMP4 
!                                                                       
      DSTRDU = -TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4+                     &
     &    NSUBI*T*(1.0-2.0*U_NUC)*(HPRIM+H/T)/AZERO)                    
!                                                                       
      DSTRDX = -TRSCAL*BRYDNS*XH*EXCLU*                                 &
     &    (MUSUBT*(DHDTDX+DHDX/T)/AZERO-                                &
     &    (T*DHDTDX+2.5*DHDX)/AZERO)                                    
!                                                                       
      DSTRDN = -TRSCAL*                                                 &
     &    (U_NUC*EXCLU*TMP4+U_NUC*EXCLU*T*(HPRIM+H/T)/AZERO)            
!                                                                       
      DSTRDT = -(BRYDNS*XH*EXCLU*((MUSUBT/T-1.5)*(HPRIM+H/T)/AZERO+     &
     &    MUSUBT*(HPPRIM+HPRIM/T-H/T**2)/AZERO-                         &
     &    (3.5*HPRIM+T*HPPRIM)/AZERO ))*TRSCAL                          
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of interior bulk !                 
!                 !      terms in the internal        !                 
!                 !      energy and entropy           !                 
!                 !      densities w.r.t. to U,X,n_i  !                 
!                 !      and T.  These are used in    !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
      S_NUC =(OVR53*UQ/T)*(TAU_NI+TAU_PI)-                              &
     &    NSUBI*((1.0-X)*ETA_NI+X*ETA_PI)                               
!                                                                       
      E_NUC = UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+      &
     &    CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM                             
!                                                                       
!                                                                       
!                    Interior particle densties                         
      NPI = X*NSUBI 
      NNI = (1.0-X)*NSUBI 
!                                                                       
      DTPIDT = 2.5*TAU_PI/T-2.25*NPI*GPI/UQ 
      DTNIDT = 2.5*TAU_NI/T-2.25*NNI*GNI/UQ 
!                                                                       
!               Derivative of interior entropy density w.r.t. T         
      DSIDT = UQ*(DTPIDT+DTNIDT)/T 
!                                                                       
!               Derivative of interior internal energy density w.r.t. T 
      DEIDT = T*DSIDT 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Derivatives of eta's w.r.t. X and NSUBI            
      DETPDX = GPI/X 
      DETNDX = -GNI/(1.0-X) 
      DETPDN = GPI/NSUBI 
      DETNDN = GNI/NSUBI 
!                                                                       
!                    Derivatives of Tau's w.r.t. X and NSUBI            
      DTPIDX = 1.5*T*NPI*DETPDX/UQ 
      DTNIDX = 1.5*T*NNI*DETNDX/UQ 
      DTPDNI = 1.5*T*NPI*DETPDN/UQ 
      DTNDNI = 1.5*T*NNI*DETNDN/UQ 
!                                                                       
!                                                                       
!                                                                       
!           Derivative of interior entropy density w.r.t. X             
      DSIDX = OVR53*UQ*(DTPIDX+DTNIDX)/T-NSUBI*(ETA_PI-ETA_NI)-         &
     &    NSUBI*((1.0-X)*DETNDX+X*DETPDX)                               
!                                                                       
!           Derivative of interior internal energy density w.r.t. X     
      DEIDX = UQ*(DTPIDX+DTNIDX)+                                       &
     &    (NSUBI**2)*4.0*BB*(1.0-2.0*X)+NSUBI*DELTAM                    
!                                                                       
!                                                                       
!           Derivative of interior entropy density w.r.t. NSUBI         
      DSIDN = OVR53*UQ*(DTPDNI+DTNDNI)/T-((1.0-X)*ETA_NI+X*ETA_PI)-     &
     &    NSUBI*((1.0-X)*DETNDN+X*DETPDN)                               
!                                                                       
!                                                                       
!           Derivative of interior internal energy density w.r.t. NSUBI 
      DEIDN = UQ*(DTPDNI+DTNDNI)+2.0*NSUBI*(AA+4.0*BB*X*(1.0-X))+       &
     &    CC*(1.0+DD)*(NSUBI**DD)+X*DELTAM                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of exterior bulk !                 
!                 !      nucleon internal energy &    !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO 
!                                                                       
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*                                &
     &((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)  
!                                                                       
!                   Derivative of exterior entropy density w.r.t. T     
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-                   &
     &     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T                            
!                                                                       
      DEODT = T*DSODT 
!                                                                       
!                    Derivatives of exterior particle densities w.r.t.  
!                    Temperature (ETA's fixed)                          
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                                                                       
      DNPDEP = NPOUT/GPO 
      DNNDEN = NNOUT/GNO 
!                                                                       
      DTPDEP = 1.5*T*NPOUT/UQ 
      DTNDEN = 1.5*T*NNOUT/UQ 
!                                                                       
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP 
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN 
!                                                                       
!                                                                       
!                    Exterior particle potentials                       
      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD ) 
      VPOUT = EIFLAG*                                                   &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)        
!                                                                       
!                                                                       
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP 
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN 
!                                                                       
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of alpha         !                 
!                 !      particle internal energy &   !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      S_ALFA = ALFDNS*(2.5-MUALFA/T) 
!                                                                       
!                                                                       
      E_ALFA = ALFDNS*(1.5*T-BALPHA) 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                                                                       
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                                                                       
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2 
!                                                                       
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS 
!                                                                       
!                                                                       
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP 
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN 
!                                                                       
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP 
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN 
!                                                                       
      DNADEP = ALFDNS*DMADEP/T 
      DNADEN = ALFDNS*DMADEN/T 
!                                                                       
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T 
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T 
!                                                                       
      DEADEP = DNADEP*(1.5*T-BALPHA) 
      DEADEN = DNADEN*(1.5*T-BALPHA) 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      S_DENS = U_NUC*S_NUC+EXCLU*EXALFA*S_OUT+EXCLU*S_ALFA+S_SC+S_TR 
!                                                                       
      E_DENS = U_NUC*E_NUC+EXCLU*EXALFA*E_OUT+EXCLU*E_ALFA+E_SC+E_TR 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !      Temperature Derivatives     !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT 
!                                                                       
!                                                                       
      DBSDT = (DU_DT*S_NUC-                                             &
     &    DU_DT*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DT*S_OUT                  &
     &    -DU_DT*S_ALFA+                                                &
     &    U_NUC*(DSIDT+DSIDX*DX_DT+DSIDN*DNI_DT)+                       &
     &    EXCLU*EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+             &
     &    EXCLU*(DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT)+                    &
     &    DSSCDT+DSSCDU*DU_DT+DSSCDX*DX_DT+DSSCDN*DNI_DT+               &
     &    DSTRDT+DSTRDU*DU_DT+DSTRDX*DX_DT+DSTRDN*DNI_DT)/BRYDNS        
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDT = T*DBSDT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+                 &
     &    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)                 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDT = BRYDNS*(DBMUDT-DBFDT) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !       Density Derivatives        !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN 
!                                                                       
!                                                                       
      DBSDN = (DU_DN*S_NUC-                                             &
     &    DU_DN*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DN*S_OUT-DU_DN*S_ALFA+    &
     &    U_NUC*(DSIDX*DX_DN+DSIDN*DNI_DN)+                             &
     &    EXCLU*EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+                   &
     &    EXCLU*(DSADEP*DEP_DN+DSADEN*DEN_DN)+                          &
     &    DSSCDU*DU_DN+DSSCDX*DX_DN+DSSCDN*DNI_DN+                      &
     &    DSTRDU*DU_DN+DSTRDX*DX_DN+DSTRDN*DNI_DN)/BRYDNS-              &
     &    S_DENS/BRYDNS**2                                              
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDN = (DU_DN*E_NUC-                                             &
     &    DU_DN*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DN*E_OUT-DU_DN*E_ALFA+    &
     &    U_NUC*(DEIDX*DX_DN+DEIDN*DNI_DN)+                             &
     &    EXCLU*EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+                   &
     &    EXCLU*(DEADEP*DEP_DN+DEADEN*DEN_DN)+                          &
     &    DESCDU*DU_DN+DESCDX*DX_DN+DESCDN*DNI_DN+                      &
     &    DETRDU*DU_DN+DETRDX*DX_DN+DETRDN*DNI_DN)/BRYDNS-              &
     &    E_DENS/BRYDNS**2                                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDN = DBUDN-T*DBSDN 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+                        &
     &    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)                        
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !         Ye Derivatives           !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY 
!                                                                       
!                                                                       
      DBSDY = (DU_DY*S_NUC-                                             &
     &    DU_DY*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DY*S_OUT-DU_DY*S_ALFA+    &
     &    U_NUC*(DSIDX*DX_DY+DSIDN*DNI_DY)+                             &
     &    EXCLU*EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+                   &
     &    EXCLU*(DSADEP*DEP_DY+DSADEN*DEN_DY)+                          &
     &    DSSCDU*DU_DY+DSSCDX*DX_DY+DSSCDN*DNI_DY+                      &
     &    DSTRDU*DU_DY+DSTRDX*DX_DY+DSTRDN*DNI_DY)/BRYDNS               
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDY = (DU_DY*E_NUC-                                             &
     &    DU_DY*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DY*E_OUT-DU_DY*E_ALFA+    &
     &    U_NUC*(DEIDX*DX_DY+DEIDN*DNI_DY)+                             &
     &    EXCLU*EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+                   &
     &    EXCLU*(DEADEP*DEP_DY+DEADEN*DEN_DY)+                          &
     &    DESCDU*DU_DY+DESCDX*DX_DY+DESCDN*DNI_DY+                      &
     &    DETRDU*DU_DY+DETRDX*DX_DY+DETRDN*DNI_DY)/BRYDNS               
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDY = DBUDY-T*DBSDY 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+                 &
     &    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN                    
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDY = BRYDNS*(DBMUDY-DBFDY) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                End of derivatives of thermodynamic variables          
!-----------------------------------------------------------------------
!                                                                       
!                  Total derivatives                                    
!                  (Baryons+Electrons+Photons)                          
!                                                                       
      DUDT = DBUDT+DEUDT+DPUDT 
      DUDN = DBUDN+DEUDN+DPUDN 
      DUDY = DBUDY+DEUDY+DPUDY 
!                                                                       
!                                                                       
      DSDT = DBSDT+DESDT+DPSDT 
      DSDN = DBSDN+DESDN+DPSDN 
      DSDY = DBSDY+DESDY+DPSDY 
!                                                                       
!                                                                       
      DPDT = DBPDT+DEPDT+DPPDT 
      DPDN = DBPDN+DEPDN+DPPDN 
      DPDY = DBPDY+DEPDY+DPPDY 
!                                                                       
!                                                                       
      DMUDT = DBMUDT+YE*DEMUDT 
      DMUDN = DBMUDN+YE*DEMUDN 
      DMUDY = DBMUDY+YE*DEMUDY 
!                                                                       
!                Calculate the adiabatic index                          
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
!                Set the value of XPREV to X for use the next           
!                time through                                           
!                                                                       
      XPREV = X 
!                                                                       
!                Save the value of the proton density to be used        
!                by the "no nuclei" scheme on the next call             
      P_PREV = NPOUT 
!                                                                       
!                                                                       
!                Return the three internal compositional variables      
      INPVAR(2) = NSUBI 
      INPVAR(3) = ETA_PO 
      INPVAR(4) = ETA_NO 
!                                                                       
!                                                                       
!                                                                       
!                Rejoice for this routine is finished!!!!!!!            
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         ALFEOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       ALFEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         8/30/90 Modified from model 4-A                      
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU                   
!                            FSWESTY@SBAST3.SUNYSB.EDU                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL ALFEOS(INPVAR,YE,BRYDNS)                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:                                                           
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                       Function type declarations                      
!                                                                       
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFO, fhalf 
!                                                                       
!                                                                       
!                                                                       
!                         Ratio of baryon density to saturation density 
      Y = BRYDNS/NSUBS 
!                                                                       
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options are not implemente
!                         in this version)                              
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                         Calc the quantum concentration of nucleons    
      NQ = 2.36D-4*T**1.5 
!                                                                       
!                         Calc the Fermi integral coefficent            
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
      KQ = ((T/UQ)**2.5)/(2.0*PI**2) 
!                                                                       
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2)) 
!                                                                       
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ) 
!                                                                       
!                                                                       
!                                                                       
!                              Set the proton density to its old value  
      NPOUT = P_PREV 
!                                                                       
      IF(BRYDNS.GT.(0.98*2.0/(YE*V_ALFA))) THEN 
        NPOUT = YE*BRYDNS 
        NNOUT = (1.0-YE)*BRYDNS 
        NOUT = BRYDNS 
!                                                                       
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
!                                                                       
        ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
        ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
        ETA_NO = FINV12(ZNO) 
!                                                                       
        ETA_PO = FINV12(ZPO) 
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        TAU_NO = KQ*F32_NO 
        TAU_PO = KQ*F32_PO 
!                                                                       
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
        MUN = MUN_O 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
        MUPROT = MUP_O 
!                                                                       
!                              Calculate diff. of chem. potentials      
        MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
!                              Calculate the alpha particle             
!                              chemical potential                       
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
        ALFDNS = 0.0 
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
      ELSE 
!                                                                       
!                              Calculate the neutron density            
        NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+         &
     &            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/                   &
     &            (2.0-BRYDNS*YE*V_ALFA)                                
!                                                                       
!                              Calculate density of outside nucleons    
        NOUT = NPOUT+NNOUT 
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
!                                                                       
        ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
        ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
        ETA_NO = FINV12(ZNO) 
!                                                                       
        ETA_PO = FINV12(ZPO) 
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        TAU_NO = KQ*F32_NO 
        TAU_PO = KQ*F32_PO 
!                                                                       
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
        MUN = MUN_O 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
        MUPROT = MUP_O 
!                                                                       
!                              Calculate diff. of chem. potentials      
        MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
!                              Calculate the alpha particle             
!                              chemical potential                       
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
!                              Calculate density of alpha particles     
!                                                                       
        IF(ABS(MUALFA/T).LT.30.0) THEN 
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
        ELSEIF((MUALFA/T).LT.-30.0) THEN 
          ALFDNS = 0.0 
        ELSE 
          ALFDNS = 8.0*NQ*DEXP(3.0D1) 
        ENDIF 
!                                                                       
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                              Calculate "non-zeroness" of baryon       
!                              conservation equation and save the       
!                              value to be used in the finite           
!                              difference approximation of DGDPRT       
        GOLD = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
        PRTOLD = NPOUT 
!                                                                       
!                              Take a small step to get derivative      
        NPOUT = NPOUT+0.001*BRYDNS 
!                                                                       
        DO 11 I=1,30,1 
!                                                                       
!                              Calculate the neutron density            
          NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+       &
     &            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/                   &
     &            (2.0-BRYDNS*YE*V_ALFA)                                
!                                                                       
!                              Calculate density of outside nucleons    
          NOUT = NPOUT+NNOUT 
!                                                                       
          VNOUT = EIFLAG*                                               &
     &      (2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)             
!                                                                       
          VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                     &
     &      CC*(1.0+DD)*NOUT**DD+DELTAM)                                
!                                                                       
!                                                                       
          ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
          ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
          ETA_NO = FINV12(ZNO) 
!                                                                       
          ETA_PO = FINV12(ZPO) 
!                                                                       
          F32_NO = F_3_2(ETA_NO) 
!                                                                       
          F32_PO = F_3_2(ETA_PO) 
!                                                                       
          TAU_NO = KQ*F32_NO 
          TAU_PO = KQ*F32_PO 
!                                                                       
          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*                           &
     &      (AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )   
!                                                                       
          MUN_O = T*ETA_NO+VNOUT 
          MUN = MUN_O 
!                                                                       
          MUP_O = T*ETA_PO+VPOUT 
          MUPROT = MUP_O 
!                                                                       
!                              Calc difference of potentials            
          MUHAT = MUN-MUPROT 
!                                                                       
!                              Calc alpha particle chemical potentials  
          MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
!                              Calc alpha particle density              
!                                                                       
          IF(ABS(MUALFA/T).LT.30.0) THEN 
            ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
          ELSEIF((MUALFA/T).LT.-30.0) THEN 
            ALFDNS = 0.0 
          ELSE 
            ALFDNS = 8.0*NQ*DEXP(3.0D1) 
          ENDIF 
!                                                                       
!                                                                       
          EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                              Calc "non-zeroness" of baryon cons. eqn. 
          G = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
!                                                                       
!                              Calculate derivative of baryon conservati
!                              equation w.r.t. proton density by finite 
!                              diference approximation                  
          DGDPRT = (G-GOLD)/(NPOUT-PRTOLD) 
!                                                                       
!                              Calculate new Newton-Raphson step        
          DPRT = G/DGDPRT 
!                                                                       
!                              Save old value of proton density & G     
          PRTOLD = NPOUT 
          GOLD = G 
!                                                                       
!                                                                       
   13     CONTINUE 
!                                                                       
!                              Potential "new" value of proton density  
          PRTNEW = NPOUT-DPRT 
!                                                                       
!                              If new proton density is less than the   
!                              baryon density and greater than zero     
!                              then update the proton density           
          IF(PRTNEW*(BRYDNS-PRTNEW).GT.0.0) THEN 
            NPOUT = NPOUT-DPRT 
!                              Else cut the step size in half and try ag
          ELSE 
            DPRT = DPRT*0.5 
            GOTO 13 
          ENDIF 
!                                                                       
!                              If step size is small enough break out of
!                              the DO 11 loop, otherwise continue       
          IF(ABS(DPRT/NPOUT).LT.10E-11) GOTO 12 
   11   CONTINUE 
!                                                                       
!      write(*,*) 'A failed to converge; switching to F' ! take out late
        SSFLAG = 0 
        GOTO 999 
!                                                                       
!                                                                       
   12   CONTINUE 
!                                                                       
      ENDIF 
!                              Set the success flag                     
      SSFLAG = 1 
!                                                                       
!                                                                       
!                              Calc outside nucleon density             
      NOUT = NNOUT+NPOUT 
!                                                                       
!                              Calc outside nucleon fraction            
      XOUT = NPOUT/NOUT 
!                                                                       
!                              Calculate particle fractions             
      XALFA = 4.0*ALFDNS/BRYDNS 
      XPROT = EXALFA*NPOUT/BRYDNS 
      XNUT = EXALFA*NNOUT/BRYDNS 
      XH = 0.0 
!                                                                       
!                              Baryons                                  
!                                                                       
      F32_NO = F_3_2(ETA_NO) 
!                                                                       
      F32_PO = F_3_2(ETA_PO) 
!                                                                       
      TAU_PO = KQ*F32_PO 
!                                                                       
      TAU_NO = KQ*F32_NO 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Calculate internal energy of outside nucleons      
      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+                        &
     &    EIFLAG*((NOUT**2)*AA+                                         &
     &   4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )/NOUT      
!                                                                       
!                                                                       
!                                Calc alfa particle internal energy     
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA) 
!                                                                       
!                                Set nuclei internal energy to zero     
      BUNUC = 0.0 
!                                Calculate total baryon internal energy 
!                                (per baryon)                           
      BU = BUOUT+BUALFA+BUNUC 
!                                                                       
!                                                                       
!                                Calc entropy of outside nucleons       
      BSOUT = (XNUT+XPROT)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-          &
     &   NNOUT*ETA_NO-NPOUT*ETA_PO )/NOUT                               
!                                                                       
!                                Calc alpha particle entropy            
      BSALFA = 0.25*XALFA*(2.5-MUALFA/T) 
!                                                                       
!                                Set nuclei entropy to zero             
      BSNUC = 0.0 
!                                                                       
!                                Calc total baryon entropy (per baryon) 
      BS = BSOUT+BSALFA+BSNUC 
!                                                                       
!                                                                       
!                                                                       
!                                Calc outside free energy               
      BFOUT = BUOUT-T*BSOUT 
!                                Calc alpha particle free energy        
      BFALFA = BUALFA-T*BSALFA 
!                                Set nuclei free energy to zero         
      BFNUC = BUNUC-T*BSNUC 
!                                Calc total baryon free energy (per nucl
      BFTOT = BFOUT+BFALFA+BFNUC 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                Calc outside pressure                  
      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                              &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)))       
!                                                                       
!                                Calc alfa particle pressure            
      BPRALF = ALFDNS*T 
!                                                                       
!                                Set nuclei pressure to zero            
      BPRNUC = 0.0 
!                                                                       
!                                Calc total baryon pressure             
      BPRESS = BPROUT+BPRALF+BPRNUC 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                           Leptons & Photons                           
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
!                                                                       
!                           Total pressure and eng/ent per baryon       
!                                                                       
      FBARY = BFTOT+FSUBE 
      PBARY = BPRESS+EPRESS 
      MUBARY = YE*MUPROT+(1.0-YE)*MUN 
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN 
!                                                                       
      FTOT = BFTOT+FSUBE+PF 
      UTOT = BU+EU+PU 
      STOT = BS+ES+PS 
      PTOT = BPRESS+EPRESS+PPRESS 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
!                Derivatives of thermodynamic variables                 
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!c      GPO = 2.0*FHALFO(ETA_PO)                                        
!c      GNO = 2.0*FHALFO(ETA_NO)                                        
!                                                                       
!                                                                       
      GPO = 2.0*FHALF(ETA_PO) 
      GNO = 2.0*FHALF(ETA_NO) 
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of exterior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and ETA's)    !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                  Derivatives of exterior potentials                   
!                  w.r.t. particle densities                            
      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVNODP = DVPODN 
      DVNODN = DVPODP 
!                                                                       
!                                                                       
!                  Derviatives of exterior chem. pot. w.r.t. ETA's      
!                  (at fixed T)                                         
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of pressure potential w.r.t.             
!                  particle densities                                   
      DV_DPO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
      DV_DNO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
!                                                                       
!                  Derivatives of pressure potential w.r.t. ETA's       
!                  (at fixed T)                                         
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
!                  Derivatives of outside pressure w.r.t. ETA's         
!                  (at fixed T)                                         
      DPODEP = NPOUT*T+DV_DEP 
      DPODEN = NNOUT*T+DV_DEN 
!                                                                       
!                  Derivatives of alpha density w.r.t. ETA's            
!                  (at fixed T)                                         
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T 
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T 
!                                                                       
!                                                                       
!                  Derivatives of particle densities w.r.t. T           
!                  (at fixed ETA's)                                     
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
!                  Derivatives of exterior chem. pot. w.r.t. T          
!                  (at fixed ETA's)                                     
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
!                  (at fixed ETA's)                                     
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
!                  (at fixed ETA's)                                     
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                  Derivative of alpha chem. pot. w.r.t. T              
!                  (at fixed ETA's)                                     
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
!                  (at fixed ETA's)                                     
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                  Derivative of alpha particle pressure w.r.t. T       
!                  (at fixed ETA's)                                     
      DPADT = ALFDNS+T*DNADT 
!                                                                       
!                                                                       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of constraint   !                  
!                 !      and equilibrium equations   !                  
!                 !      with respect to the five    !                  
!                 !      compositional variables     !                  
!                 !      (U,x,n_i,eta_po,eta_no)     !                  
!                 !      and the three independent   !                  
!                 !      variables                   !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 1 (Baryon conservation)                       
!                                                                       
!                                                                       
!                                                                       
      DG1DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NOUT-4.0)*DNADEP 
!                                                                       
      DG1DO2 = -EXALFA*NNOUT/GNO+(V_ALFA*NOUT-4.0)*DNADEN 
!                                                                       
!                                                                       
      DG1DL1 = 1.0 
!                                                                       
      DG1DL2 = -EXALFA*(DNNODT+DNPODT)+(V_ALFA*NOUT-4.0)*DNADT 
!                                                                       
      DG1DL3 = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 2 (Charge conservation)                       
!                                                                       
!                                                                       
      DG2DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NPOUT-2.0)*DNADEP 
!                                                                       
      DG2DO2 = (V_ALFA*NPOUT-2.0)*DNADEN 
!                                                                       
!                                                                       
      DG2DL1 = YE 
!                                                                       
      DG2DL2 = -EXALFA*DNPODT+(V_ALFA*NPOUT-2.0)*DNADT 
!                                                                       
      DG2DL3 = BRYDNS 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      DET_GT = DG1DO1*DG2DO2-DG1DO2*DG2DO1 
!                                                                       
!                                                                       
      DEP_DN = (DG1DO2*DG2DL1-DG2DO2*DG1DL1)/DET_GT 
      DEN_DN = (DG2DO1*DG1DL1-DG1DO1*DG2DL1)/DET_GT 
!                                                                       
!                                                                       
      DEP_DT = (DG1DO2*DG2DL2-DG2DO2*DG1DL2)/DET_GT 
      DEN_DT = (DG2DO1*DG1DL2-DG1DO1*DG2DL2)/DET_GT 
!                                                                       
!                                                                       
      DEP_DY = (DG1DO2*DG2DL3-DG2DO2*DG1DL3)/DET_GT 
      DEN_DY = (DG2DO1*DG1DL3-DG1DO1*DG2DL3)/DET_GT 
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of exterior bulk !                 
!                 !      nucleon internal energy &    !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO 
!                                                                       
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*                                &
     &((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)  
!                                                                       
!                   Derivative of exterior entropy density w.r.t. T     
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-                   &
     &     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T                            
!                                                                       
      DEODT = T*DSODT 
!                                                                       
!                    Derivatives of exterior particle densities w.r.t.  
!                    Temperature (ETA's fixed)                          
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                                                                       
      DNPDEP = NPOUT/GPO 
      DNNDEN = NNOUT/GNO 
!                                                                       
      DTPDEP = 1.5*T*NPOUT/UQ 
      DTNDEN = 1.5*T*NNOUT/UQ 
!                                                                       
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP 
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN 
!                                                                       
!                                                                       
!                    Exterior particle potentials                       
      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD ) 
      VPOUT = EIFLAG*                                                   &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)        
!                                                                       
!                                                                       
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP 
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN 
!                                                                       
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of alpha         !                 
!                 !      particle internal energy &   !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      S_ALFA = ALFDNS*(2.5-MUALFA/T) 
!                                                                       
!                                                                       
      E_ALFA = ALFDNS*(1.5*T-BALPHA) 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                                                                       
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                                                                       
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2 
!                                                                       
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS 
!                                                                       
!                                                                       
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP 
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN 
!                                                                       
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP 
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN 
!                                                                       
      DNADEP = ALFDNS*DMADEP/T 
      DNADEN = ALFDNS*DMADEN/T 
!                                                                       
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T 
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T 
!                                                                       
      DEADEP = DNADEP*(1.5*T-BALPHA) 
      DEADEN = DNADEN*(1.5*T-BALPHA) 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      S_DENS = EXALFA*S_OUT+S_ALFA 
!                                                                       
      E_DENS = EXALFA*E_OUT+E_ALFA 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !      Temperature Derivatives     !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT 
!                                                                       
!                                                                       
      DBSDT = (-V_ALFA*DNA_DT*S_OUT+                                    &
     &    EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+                   &
     &    (DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT) )/BRYDNS                  
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDT = T*DBSDT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+                 &
     &    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)                 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDT = BRYDNS*(DBMUDT-DBFDT) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !       Density Derivatives        !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN 
!                                                                       
!                                                                       
      DBSDN = (-V_ALFA*DNA_DN*S_OUT+                                    &
     &    EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+                         &
     &   (DSADEP*DEP_DN+DSADEN*DEN_DN) )/BRYDNS-S_DENS/BRYDNS**2        
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBUDN = (-V_ALFA*DNA_DN*E_OUT+                                    &
     &    EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+                         &
     &   (DEADEP*DEP_DN+DEADEN*DEN_DN) )/BRYDNS-E_DENS/BRYDNS**2        
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDN = DBUDN-T*DBSDN 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+                        &
     &    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)                        
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !         Ye Derivatives           !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY 
!                                                                       
!                                                                       
      DBSDY = (-V_ALFA*DNA_DY*S_OUT+                                    &
     &    EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+                         &
     &   (DSADEP*DEP_DY+DSADEN*DEN_DY) )/BRYDNS                         
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBUDY = (-V_ALFA*DNA_DY*E_OUT+                                    &
     &    EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+                         &
     &   (DEADEP*DEP_DY+DEADEN*DEN_DY) )/BRYDNS                         
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDY = DBUDY-T*DBSDY 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+                 &
     &    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN                    
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDY = BRYDNS*(DBMUDY-DBFDY) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                End of derivatives of thermodynamic variables          
!-----------------------------------------------------------------------
!                                                                       
!                  Total derivatives                                    
!                  (Baryons+Electrons+Photons)                          
!                                                                       
      DUDT = DBUDT+DEUDT+DPUDT 
      DUDN = DBUDN+DEUDN+DPUDN 
      DUDY = DBUDY+DEUDY+DPUDY 
!                                                                       
!                                                                       
      DSDT = DBSDT+DESDT+DPSDT 
      DSDN = DBSDN+DESDN+DPSDN 
      DSDY = DBSDY+DESDY+DPSDY 
!                                                                       
!                                                                       
      DPDT = DBPDT+DEPDT+DPPDT 
      DPDN = DBPDN+DEPDN+DPPDN 
      DPDY = DBPDY+DEPDY+DPPDY 
!                                                                       
!                                                                       
      DMUDT = DBMUDT+YE*DEMUDT 
      DMUDN = DBMUDN+YE*DEMUDN 
      DMUDY = DBMUDY+YE*DEMUDY 
!                                                                       
!                                                                       
!                  Adiabatic index                                      
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
      INPVAR(2) = NSUBS 
      INPVAR(3) = ETA_PO 
      INPVAR(4) = ETA_NO 
!                                                                       
!                                                                       
!                           Approximate the nuclear density             
      NSUBI = NSUBS 
!                                                                       
!                           Use 0.45 as the nuclear proton fraction     
      X = 0.45 
!                                                                       
!                           Save the proton number density for use      
!                           as the initial guess on next call           
      P_PREV = NPOUT 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         MAXWEL.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MAXWEL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/13/90                                              
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU                   
!                                                                       
!                                                                       
!    CALL LINE:    CALL MAXWEL(INPVAR,YE,BRYDNS)                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:                                                           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION OUTVAR(4) 
      DOUBLE PRECISION DPN_DT, DPN_DN, DPN_DY 
      DOUBLE PRECISION DSN_DT, DSN_DN, DSN_DY 
      DOUBLE PRECISION DSB_DT, DSB_DN, DSB_DY 
      DOUBLE PRECISION DMU_DT, DMU_DN, DMU_DY 
      DOUBLE PRECISION DPHADT, DPHADY, DELDNS 
      DOUBLE PRECISION N_XH, N_XA, N_XN, N_XP, B_XA, B_XN, B_XP 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
!                                                                       
!                   Set the temperature                                 
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                   Calculate and save chem. pot. and thermodynamic     
!                   quantaties from low end of two phase region         
      CALL NUCEOS(INPVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
      IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
        CALL RESET(INPVAR,YE,LOWDNS,OUTVAR) 
        OUTVAR(1) = INPVAR(1) 
        CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
        IF(SSFLAG.NE.1) THEN 
          OUTVAR(2) = 0.155 
          OUTVAR(3) = -15.0 
          OUTVAR(4) = -20.0 
          CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
        ELSE 
          INPVAR(2) = OUTVAR(2) 
          INPVAR(3) = OUTVAR(3) 
          INPVAR(4) = OUTVAR(4) 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      PRLOW = PTOT-PPRESS 
      S_LOW = STOT-PS 
      F_LOW = FTOT-PF 
      MUTLOW = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE) 
      MUELOW = MUSUBE 
      MUHLOW = MUHAT 
!                                                                       
      DPN_DT = DPDT 
      DPN_DN = DPDN 
      DPN_DY = DPDY 
!                                                                       
      DMU_DT = DMUDT 
      DMU_DN = DMUDN 
      DMU_DY = DMUDY 
!                                                                       
      DSN_DT = DSDT-DPSDT 
      DSN_DN = DSDN 
      DSN_DY = DSDY 
!                                                                       
      N_XH = XH 
      N_XA = XALFA 
      N_XP = XPROT 
      N_XN = XNUT 
!                                                                       
!                                                                       
      IF(SSFLAG.NE.1) THEN 
        WRITE(*,*) 'MAXWEL:  Nuclear EOS failed at try:' 
        WRITE(*,*) T,LOWDNS,YE 
        WRITE(*,*) INPVAR 
        GOTO 999 
      ENDIF 
!                   Calculate and save chem. pot. and thermodynamic     
!                   quantaties from high end of two phase region        
      CALL ALFEOS(INPVAR,YE,HIDNS,P_PREV,SSFLAG) 
!                                                                       
      PRHI = PTOT-PPRESS 
      S_HI = STOT-PS 
      F_HI = FTOT-PF 
      MUTHI = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE) 
      MUEHI = MUSUBE 
      MUHHI = MUHAT 
!                                                                       
!                                                                       
      DSB_DT = DSDT-DPSDT 
      DSB_DN = DSDN 
      DSB_DY = DSDY 
!                                                                       
!                                                                       
      B_XA = XALFA 
      B_XP = XPROT 
      B_XN = XNUT 
!                                                                       
!                                                                       
      IF(SSFLAG.NE.1) THEN 
        WRITE(*,*) 'MAXWEL:  Alfa EOS failed at try:' 
        WRITE(*,*) T,HIDNS,YE 
        WRITE(*,*) INPVAR 
        GOTO 999 
      ENDIF 
!                                                                       
!                   Calculate "average" chem. pot. and pressure         
!                   in order to avoid numerical problems                
      MUTILD = (MUTLOW+MUTHI)/2.0 
      PRTILD = (PRLOW+PRHI)/2.0 
!                                                                       
!                   Calculate phase fraction                            
      PHASEF = (BRYDNS-LOWDNS)/(HIDNS-LOWDNS) 
!                                                                       
!                                                                       
!                   Electron number density                             
      NSUBE = BRYDNS*YE 
!                                                                       
!                   Call electron EOS to determine the                  
!                   electron chemical potential                         
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
      MUHAT = MUSUBE+(1.0-PHASEF)*(MUHLOW-MUELOW)+PHASEF*(MUHHI-MUEHI) 
!                                                                       
      MUN = MUTILD+YE*(MUHAT-MUSUBE) 
!                                                                       
      MUPROT = MUN-MUHAT 
!                                                                       
!                   Calculate thermodynamic quantities                  
!                                                                       
      STOT = ((1.0-PHASEF)*S_LOW*LOWDNS+PHASEF*S_HI*HIDNS)/BRYDNS+PS 
!                                                                       
      FTOT = (LOWDNS*F_LOW+MUTILD*(BRYDNS-LOWDNS))/BRYDNS+PF 
!                                                                       
      UTOT = FTOT+T*STOT+PU 
!                                                                       
      PTOT = PRTILD+PPRESS 
!                                                                       
!                                                                       
      XH = (1.0-PHASEF)*N_XH 
      XALFA = (1.0-PHASEF)*N_XA 
      XNUT = (1.0-PHASEF)*N_XN 
      XPROT = (1.0-PHASEF)*N_XP 
      XALFA2 = PHASEF*B_XA 
      XNUT2 = PHASEF*B_XN 
      XPROT2 = PHASEF*B_XP 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DELDNS = HIDNS-LOWDNS 
!                                                                       
!                                                                       
      DPHADT = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DT-           &
     &    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DT                            
!                                                                       
      DPDT = DPN_DT+DPN_DN*DNL_DT 
      DMUDT = DMU_DT+DMU_DN*DNL_DT 
      DSDT = (1.0-PHASEF)*LOWDNS*(DSN_DT+DSN_DN*DNL_DT)/BRYDNS+         &
     & (1.0-PHASEF)*S_LOW*DNL_DT/BRYDNS-LOWDNS*S_LOW*DPHADT/BRYDNS+     &
     &    (DPHADT*S_HI*HIDNS+PHASEF*DNH_DT*S_HI+                        &
     &    PHASEF*HIDNS*(DSB_DT+DSB_DN*DNH_DT))/BRYDNS+DPSDT             
      DUDT = DMUDT-DPDT/BRYDNS+STOT+T*DSDT 
!                                                                       
!                                                                       
      DPDN = 0.0 
      DMUDN = 0.0 
      DSDN = -DPDT/BRYDNS**2 
      DUDN = (LOWDNS*(MUTILD-FTOT)/BRYDNS**2)+T*DSDN 
!                                                                       
!                                                                       
      DPHADY = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DY-           &
     &    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DY                            
!                                                                       
      DPDY = DPN_DY+DPN_DN*DNL_DY 
      DMUDY = DMU_DY+DMU_DN*DNL_DY 
      DSDY = (1.0-PHASEF)*LOWDNS*(DSN_DY+DSN_DN*DNL_DY)/BRYDNS+         &
     & (1.0-PHASEF)*S_LOW*DNL_DY/BRYDNS-LOWDNS*S_LOW*DPHADY/BRYDNS+     &
     &    (DPHADY*S_HI*HIDNS+PHASEF*DNH_DY*S_HI+                        &
     &    PHASEF*HIDNS*(DSB_DY+DSB_DN*DNH_DY))/BRYDNS                   
      DUDY = DMUDY-DPDY/BRYDNS+T*DSDY 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!             Adiabatic index                                           
!             (Note that the first term vanishes in this expression)    
      GAM_S = T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EOSLOG.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EOSLOG                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/15/90                                             
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)                 
!                                                                       
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!                                                                       
!                                                                       
!    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B           
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS (includes alpha's)           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC, MAXWEL.INC                            
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EOSLOG(INPVAR,YE,BRYDNS,EOSFLG) 
!                                                                       
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                       This include file contains all variable         
!                       declarations.  NOTE:: no implicit typing        
!                       scheme is followed in this code; if you         
!                       have any doubt as to a variables type CHECK     
!                       IT!!!!.  Also note that ALL variables are       
!                       declared explicitly.                            
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
      DOUBLE PRECISION NLOW, NHI, N_CUT, TEMP_1, TEMP_2, T_BNDY 
!                                                                       
      DOUBLE PRECISION LMM, LMP, LPM, LPP 
      DOUBLE PRECISION DNDY1, DNDY2 
!                                                                       
!                                                                       
!                                                                       
   10 CONTINUE 
!                                                                       
!                                                                       
!                         Set T equal to the input variable (any calls  
!                         with entropy or internal energy should go     
!                         the the EOS_M4B subroutine)                   
!                                                                       
      T = INPVAR(1) 
!                                                                       
!-----------------------------------------------------------------------
!         code to figure out the boundaries from the tables             
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF(YE.GT.Y_HI) THEN 
!                         Ye is too large for EOS                       
!                                                                       
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time' 
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_HI,' instead' 
        YE = Y_HI 
        GOTO 10 
!                                                                       
      ELSEIF(YE.GT.Y_LOW) THEN 
!                         Calculate high and low boundary densities     
!                         for the Maxwell construction                  
!                                                                       
!----------------------------------------------------------             
!           Calc Ye index                                               
!----------------------------------------------------------             
!                                                                       
        YFRAC = (YE-Y_LOW)/(Y_HI-Y_LOW) 
        J_MXWL = INT(YFRAC*(NUMYE-1))+1 
        DELT_Y = (Y_HI-Y_LOW)/DBLE(NUMYE-1) 
!                                                                       
        YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
        YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
!                                                                       
!                                                                       
        IF((YE.GE.YMINUS).AND.(YE.LE.YPLUS)) THEN 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
        ELSEIF(YE.GT.YPLUS) THEN 
          J_MXWL = J_MXWL+1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ELSE 
          J_MXWL = J_MXWL-1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ENDIF 
!                                                                       
!                                                                       
        IF(J_MXWL.GT.(NUMYE-1)) THEN 
          J_MXWL = NUMYE-1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ENDIF 
!                                                                       
!                                                                       
        YINTRP = (YE-YMINUS)/(YPLUS-YMINUS) 
!                                                                       
!----------------------------------------------------------             
!           Calc T index                                                
!----------------------------------------------------------             
!                                                                       
!                                                                       
        TFRAC = (T-T_LOW)/(T_HI-T_LOW) 
        I_MXWL = INT(TFRAC*(NUMTMP-1))+1 
        DELT_T = (T_HI-T_LOW)/DBLE(NUMTMP-1) 
!                                                                       
        TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
        TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
!                                                                       
!                                                                       
        IF((T.GT.TMINUS).AND.(T.LE.TPLUS)) THEN 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ELSEIF(T.GT.TPLUS) THEN 
          I_MXWL = I_MXWL+1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ELSE 
          I_MXWL = I_MXWL-1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ENDIF 
!                                                                       
!                                                                       
        IF(I_MXWL.GT.(NUMTMP-1)) THEN 
          I_MXWL = NUMTMP-1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ENDIF 
!                                                                       
!                                                                       
        TINTRP = (T-TMINUS)/(TPLUS-TMINUS) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                Find the temperature and density at the top of the     
!                Maxwel construction                                    
!                                                                       
!C      T_MXWL = YINTRP*(T_H(J_MXWL+1)-T_H(J_MXWL))+T_H(J_MXWL)         
!C      D_MXWL = YINTRP*(D_H(J_MXWL+1)-D_H(J_MXWL))+D_H(J_MXWL)         
        T_MXWL = DMIN1(T_H(J_MXWL+1),T_H(J_MXWL)) 
        IF(T_H(J_MXWL+1).GT.T_H(J_MXWL)) THEN 
          D_MXWL = D_H(J_MXWL) 
        ELSE 
          D_MXWL = D_H(J_MXWL+1) 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
!--------------------------------------------------------------------   
!            Interpolate to get Maxwell construction densities          
!--------------------------------------------------------------------   
!                                                                       
!                                                                       
!                                                                       
        DNS_1 = YINTRP*(BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))+ &
     &               BRYLOW(I_MXWL,J_MXWL)                              
        DNS_2 = YINTRP*                                                 &
     &        (BRYLOW(I_MXWL+1,J_MXWL+1)-BRYLOW(I_MXWL+1,J_MXWL))+      &
     &               BRYLOW(I_MXWL+1,J_MXWL)                            
!                                                                       
        LOWDNS = TINTRP*(DNS_2-DNS_1)+DNS_1 
!                                                                       
!                Derivative of lower density w.r.t. T                   
        DNL_DT = (DNS_2-DNS_1)/DELT_T 
!                                                                       
        DNDY1 = (BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))/DELT_Y 
        DNDY2 = (BRYLOW(I_MXWL+1,J_MXWL+1)-                             &
     &      BRYLOW(I_MXWL+1,J_MXWL))/DELT_Y                             
        DNL_DY = TINTRP*(DNDY2-DNDY1)+DNDY1 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        IF(YE.GT.Y_CUT) THEN 
!                                                                       
          DNS_1 = YINTRP*                                               &
     &        (BRYHI(I_MXWL,J_MXWL+1)-BRYHI(I_MXWL,J_MXWL))+            &
     &        BRYHI(I_MXWL,J_MXWL)                                      
          DNS_2 = YINTRP*                                               &
     &        (BRYHI(I_MXWL+1,J_MXWL+1)-BRYHI(I_MXWL+1,J_MXWL))+        &
     &               BRYHI(I_MXWL+1,J_MXWL)                             
!                                                                       
          HIDNS = TINTRP*(DNS_2-DNS_1)+DNS_1 
!                                                                       
!                Derivative of higher density w.r.t. T                  
          DNH_DT = (DNS_2-DNS_1)/DELT_T 
!                                                                       
!                                                                       
        DNDY1 = (BRYHI(I_MXWL,J_MXWL+1)-                                &
     &      BRYHI(I_MXWL,J_MXWL))/DELT_Y                                
        DNDY2 = (BRYHI(I_MXWL+1,J_MXWL+1)-                              &
     &      BRYHI(I_MXWL+1,J_MXWL))/DELT_Y                              
        DNH_DY = TINTRP*(DNDY2-DNDY1)+DNDY1 
!                                                                       
!                                                                       
        ELSE 
          HIDNS = LOWDNS 
        ENDIF 
!                                                                       
!                                                                       
!--------------------------------------------------------------------   
!--------------------------------------------------------------------   
!                                                                       
!                       Ye is too low                                   
      ELSE 
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time' 
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_LOW,' instead' 
        YE = Y_LOW 
        GOTO 10 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DLTLN1 = (LNCUT-LNLOW)/DBLE(NUMLOW-1) 
      DLTLN2 = (LNHI-LNCUT)/DBLE(NUMHI-1) 
!                                                                       
!                                                                       
      NLOW = 10.0**LNLOW 
      NHI = 10.0**LNHI 
      N_CUT = 10.0**LNCUT 
      LOGBRY = DLOG10(BRYDNS) 
      LOGBCH = LOGBRY 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!----------------------------------------------------------             
!           Calc T index                                                
!----------------------------------------------------------             
!                                                                       
!                                                                       
      IF(LOGBRY.GE.LNHI) THEN 
        I_BD = NBPNTS 
        I_BNDY = NBPNTS 
        T_BNDY = YINTRP*                                                &
     &           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+       &
     &            LBOUND(I_BNDY,J_BNDY)                                 
        GOTO 70 
      ELSEIF((LOGBRY.LT.LNHI).AND.(LOGBRY.GT.LNCUT)) THEN 
!                                                                       
        I_BD = INT((LOGBRY-LNCUT)/DLTLN2)+NUMLOW 
        LNMINS = LNCUT+DBLE(I_BD-NUMLOW)*DLTLN2 
        LNPLUS = LNCUT+DBLE(I_BD-NUMLOW+1)*DLTLN2 
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN 
          I_BNDY = I_BD 
        ELSEIF(LOGBCH.GT.LNPLUS) THEN 
          I_BD = I_BD+1 
          I_BNDY = I_BD 
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
        ELSE 
          I_BD = I_BD-1 
          I_BNDY = I_BD 
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
        ENDIF 
!                                                                       
      ELSEIF((LOGBRY.LE.LNCUT).AND.(LOGBRY.GT.LNLOW)) THEN 
!                                                                       
        I_BD = INT((LOGBRY-LNLOW)/DLTLN1)+1 
        LNMINS = LNLOW+DBLE(I_BD-1)*DLTLN1 
        LNPLUS = LNLOW+DBLE(I_BD)*DLTLN1 
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN 
          I_BNDY = I_BD 
        ELSEIF(LOGBCH.GT.LNPLUS) THEN 
          I_BD = I_BD+1 
          I_BNDY = I_BD 
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1 
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1 
        ELSE 
          I_BD = I_BD-1 
          I_BNDY = I_BD 
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1 
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      IF(I_BNDY.GT.(NBPNTS-1)) THEN 
        I_BD = NBPNTS-1 
        I_BNDY = I_BD 
        LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
        LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
      LMM = LBOUND(I_BNDY,J_BNDY) 
      LPM = LBOUND(I_BNDY+1,J_BNDY) 
      LMP = LBOUND(I_BNDY,J_BNDY+1) 
      LPP = LBOUND(I_BNDY+1,J_BNDY+1) 
!                                                                       
      LNFRAC = (LOGBCH-LNMINS)/(LNPLUS-LNMINS) 
!                                                                       
!                Interpolate in Ye first                                
!                                                                       
      TEMP_1 = YINTRP*                                                  &
     &           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+       &
     &            LBOUND(I_BNDY,J_BNDY)                                 
      TEMP_2 = YINTRP*                                                  &
     &        (LBOUND(I_BNDY+1,J_BNDY+1)-LBOUND(I_BNDY+1,J_BNDY))+      &
     &               LBOUND(I_BNDY+1,J_BNDY)                            
!                                                                       
!                Interpolate in density between the two Ye              
!                interpolated values                                    
!                                                                       
      T_BNDY = LNFRAC*(TEMP_2-TEMP_1)+TEMP_1 
!                                                                       
!                                                                       
!----------------------------------------------------------             
!----------------------------------------------------------             
!                                                                       
   70 CONTINUE 
!                                                                       
      TCHK_B = 1.01*T_BNDY 
      TCHK_N = 0.95*T_BNDY 
!                                                                       
      IF((LMM.GE.LPM).OR.(LMP.GT.LPP)) THEN 
        TCHK_N = DMAX1(0.0D0,DMIN1(0.95*TCHK_N,T_BNDY-3.0)) 
      ENDIF 
!                                                                       
!-----------------------------------------------------------------------
!               EOS Logic                                               
!-----------------------------------------------------------------------
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                     If T is below the maximum for maxwel construction 
      IF(T.LT.T_MXWL) THEN 
!                       If rho is greater than the upper max. con.      
!                       density the use the bulk EOS                    
        IF(BRYDNS.GT.HIDNS) THEN 
          EOSFLG = 3 
!                       Else if rho is greater than the lower max. con. 
!                       density then                                    
        ELSEIF(BRYDNS.GT.LOWDNS) THEN 
!                         If Ye is large enough to have a signifigant   
!                         max con then use the maxwell con. EOS         
          IF(YE.GT.Y_CUT) THEN 
            EOSFLG = 4 
!                         Otherwise use the bulk EOS                    
          ELSE 
            EOSFLG = 3 
          ENDIF 
!                                                                       
!                       If density is greater than the minimum          
!                       Maxwell con. density, then we know that we are  
!                       in the Nuclear EOS density                      
        ELSEIF(BRYDNS.GT.D_MXWL) THEN 
          EOSFLG = 2 
!                                                                       
!                                                                       
!                       Otherwise check the Boundary table              
        ELSE 
!                                                                       
!                         If T is well below the phase boundary curve   
!                         then use the nuclear EOS                      
          IF(T.LT.TCHK_N) THEN 
            EOSFLG = 2 
!                         Otherwise if T is near the boundary, first    
!                         try the nuclear EOS and if not successfull    
!                         then use the bulk EOS                         
          ELSEIF(T.LT.TCHK_B) THEN 
            EOSFLG = 1 
          ELSE 
!                         Otherwise T is well above the boundary so     
!                         use the bulk EOS                              
            EOSFLG = 3 
          ENDIF 
        ENDIF 
!                                                                       
!                     Otherwise T is above the maximum for a maxwell    
!                     construction                                      
      ELSE 
!                       If density is greater than that at the top of   
!                       the maxwell construction then use the bulk EOS  
        IF(BRYDNS.GT.D_MXWL) THEN 
          EOSFLG = 3 
!                                                                       
!                       Otherwise density is below the maxwell con.     
        ELSE 
!                                                                       
!                         If T is well below the phase boundary curve   
!                         then use the nuclear EOS                      
          IF(T.LT.TCHK_N) THEN 
            EOSFLG = 2 
!                                                                       
!                         Otherwise if T is near the phase boundary     
!                         curve then try the nuclear EOS and if not     
!                         successfull then use the bulk EOS             
          ELSEIF(T.LT.TCHK_B) THEN 
            EOSFLG = 1 
!                                                                       
!                         Otherwise T is well above the phase boundary  
!                         curve so use the bulk EOS                     
          ELSE 
            EOSFLG = 3 
          ENDIF 
        ENDIF 
      ENDIF 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                         Done with EOS logic so return EOSFLG          
!-----------------------------------------------------------------------
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         RESET.FOR                                            
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       RESET                                                
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/21/90                                             
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)                  
!                                                                       
!                                                                       
!    INPUTS:       INPVAR = TEMP, NSUBI, ETA_PO, ETA_NO                 
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!                                                                       
!                                                                       
!    OUTPUTS:      OUTVAR = ARRAY OF LENGTH 4 CONTAINING RESET VALUES   
!                  FOR THE INITIAL GUESSES                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES: NONE                                                
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE RESET(INPVAR,YE,BRYDNS,OUTVAR) 
!                                                                       
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
!                      Subroutine parameters                            
!                                                                       
      DOUBLE PRECISION INPVAR(4), OUTVAR(4), YE, BRYDNS 
!                                                                       
!                                                                       
!                      Local variables                                  
!                                                                       
      DOUBLE PRECISION ZPG, ZNG, ETA_PG, ETA_NG, PI, UQ, MQ, T, EFRAC 
!                                                                       
!                      Functions                                        
!                                                                       
      DOUBLE PRECISION FINV12 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      T = INPVAR(1) 
!                                                                       
!                                                                       
      PI = 3.1415927 
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
!                                                                       
      EFRAC = 0.5*YE 
!                                                                       
      ZNG = 2.0*(PI**2)*BRYDNS*(1.0-EFRAC)/MQ 
!                                                                       
      ZPG = 2.0*(PI**2)*BRYDNS*EFRAC/MQ 
!                                                                       
      ETA_NG = FINV12(ZNG) 
!                                                                       
      ETA_PG = FINV12(ZPG) 
!                                                                       
      OUTVAR(1) = INPVAR(1) 
      OUTVAR(2) = INPVAR(2) 
      OUTVAR(3) = ETA_PG 
      OUTVAR(4) = ETA_NG 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         ALOADMX.FOR                                          
!    MODULE:       LOADMX                                               
!    TYPE:         LOADMX                                               
!                                                                       
!    PURPOSE:      LOAD THE LOOK-UP TABLE FOR THE MAXWELL CONSTRUCTION  
!                                                                       
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/16/90                                              
!                                                                       
!    CALL LINE:    CALL LOADMX                                          
!                                                                       
!    INPUTS:       N/A                                                  
!                                                                       
!    OUTPUTS       N/A                                                  
!                                                                       
!    SUBROUTINE CALLS: EOS_M4A                                          
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC, MAXWEL.INC                            
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE LOADMX() 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
      INTEGER NTMP, NYE, NYE2, NUM_BP 
      INTEGER LUN1, LUN2, KK, KMIN 
      PARAMETER(LUN1=44,LUN2=45) 
!                                                                       
!                                                                       
      INTEGER FNML1, FNML2 
      CHARACTER*60 FNAME1, FNAME2 
!                                                                       
      DOUBLE PRECISION N_SM, SYMM_M, COMP_M, BINDEM, SYM_SM, SIG_SM 
      DOUBLE PRECISION N_SB, SYMM_B, COMP_B, BINDEB, SYM_SB, SIG_SB 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!******below commented out by MBD, always read in max180.atb            
!******and bd180.atb                                                    
!      CALL GETFNM(FNAME1,FNML1,'Enter ASCII Maxwell fname:',26)        
!                                                                       
!      CALL GETFNM(FNAME2,FNML2,'Enter ASCII boundary fname:',27)       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!        Read the file Maxwell construction data file                   
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
      OPEN(UNIT=LUN1,FILE='max180.atb',STATUS='OLD') 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN1,*) N_SM, SYMM_M 
      READ(LUN1,*) COMP_M,BINDEM 
      READ(LUN1,*) SYM_SM, SIG_SM 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN1,*) NTMP,NYE 
      READ(LUN1,*) T_LOW,T_HI 
      READ(LUN1,*) Y_LOW,Y_HI 
!                                                                       
!                                                                       
!                                                                       
      IF((NTMP.NE.NUMTMP).OR.(NYE.NE.NUMYE)) THEN 
        WRITE(*,*) 'LOADMX:  MXWL TABLE IS INCOMPATIBLE WITH ARRAYS' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      DO 101 J=1,NUMYE,1 
        DO 100 I=1,NUMTMP,3 
          KMIN = MIN0(I+2,NUMTMP) 
          READ(LUN1,*) (BRYLOW(KK,J),KK=I,KMIN,1) 
  100   CONTINUE 
  101 END DO 
!                                                                       
!                                                                       
      DO 103 J=1,NUMYE,1 
        DO 102 I=1,NUMTMP,3 
          KMIN = MIN0(I+2,NUMTMP) 
          READ(LUN1,*) (BRYHI(KK,J),KK=I,KMIN,1) 
  102   CONTINUE 
  103 END DO 
!                                                                       
!                                                                       
!                                                                       
      DO 104 I=1,NUMYE,3 
        KMIN = MIN0(I+2,NUMYE) 
        READ(LUN1,*) (T_H(KK),KK=I,KMIN,1) 
  104 END DO 
!                                                                       
!                                                                       
      DO 105 I=1,NUMYE,3 
        KMIN = MIN0(I+2,NUMYE) 
        READ(LUN1,*) (D_H(KK),KK=I,KMIN,1) 
  105 END DO 
!                                                                       
      READ(LUN1,*) YCUT 
!                                                                       
!                                                                       
      CLOSE(UNIT=LUN1,STATUS='KEEP') 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  MAXWELL CON. TABLE IS INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!        Read the file Boundary data file                               
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
      OPEN(UNIT=LUN2,FILE='bd180.atb',STATUS='OLD') 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN2,*) N_SB,SYMM_B 
      READ(LUN2,*) COMP_B,BINDEB 
      READ(LUN2,*) SYM_SB,SIG_SB 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN2,*) NUM_BP,NYE2 
      READ(LUN2,*) LNL,LNH,LNC 
      READ(LUN2,*) Y_LOW2,Y_HI2 
!                                                                       
!                                                                       
      IF((NBPNTS.NE.NUM_BP).OR.(NYE2.NE.NUMYE)) THEN 
        WRITE(*,*) 'LOADMX:  BNDY TABLE IS INCOMPATIBLE WITH ARRAYS' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(LNL-LNLOW).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  LOWER END OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      IF(ABS(LNH-LNHI).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  UPPER END OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      IF(ABS(LNC-LNCUT).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  MID CUT OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(Y_LOW-Y_LOW2).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  LOWER YE LIMITS ARE INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(Y_HI-Y_HI2).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  UPPER YE LIMITS ARE INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      DO 201 J=1,NUMYE,1 
        DO 200 I=1,NBPNTS,3 
          KMIN = MIN0(I+2,NBPNTS) 
          READ(LUN2,*) (LBOUND(KK,J),KK=I,KMIN,1) 
  200   CONTINUE 
  201 END DO 
!                                                                       
!                                                                       
      DO 203 J=1,NUMYE,1 
        DO 202 I=1,NBPNTS,3 
          KMIN = MIN0(I+2,NBPNTS) 
          READ(LUN2,*) (UBOUND(KK,J),KK=I,KMIN,1) 
  202   CONTINUE 
  203 END DO 
!                                                                       
!                                                                       
!                                                                       
      CLOSE(UNIT=LUN2,STATUS='KEEP') 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  BOUNDARY TABLE IS INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                  All arrays are now loaded so return                  
!-----------------------------------------------------------------------
!                                                                       
      N_S = N_SM 
      NSUBS = N_SM 
      SYMM = SYMM_M 
      COMP = COMP_M 
      BIND_E = BINDEM 
      SYM_S = SYM_SM 
      SIG_S = SIG_SM 
!                                                                       
      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23 
      DD = (COMP+2.0*SKYRMC)/(3.0*SKYRMC+9.0*BIND_E) 
      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S 
      AA = (OVR23*SKYRMC-DD*(SKYRMC+BIND_E))/(N_S*(DD-1.0))-BB 
      CC = (COMP+2.0*SKYRMC)/(9.0*DD*(DD-1.0)*N_S**DD) 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  SKYRME PARAMETERS FOR THIS RUN ARE:>>' 
      WRITE(*,*) 'ABCD: ',AA,BB,CC,DD 
      WRITE(*,*) ' Satur. density, symmetry engy, & compression mod.:' 
      WRITE(*,*) N_SM, SYMM_M, COMP_M 
      WRITE(*,*) N_SB, SYMM_B, COMP_B 
      WRITE(*,*) ' Binding engy, surf. symm. engy, & surface tension:' 
      WRITE(*,*) BINDEM,SYM_SM,SIG_SM 
      WRITE(*,*) BINDEB,SYM_SB,SIG_SB 
!                                                                       
      WRITE(*,*) 
!                                                                       
!                                                                       
      CALL INITFERM() 
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX: FERMI INTEGRAL TABLES ARE INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       GETFNM.FOR                                           
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         8/5/89                                               
!                                                                       
!    PURPOSE:      OBTAINS A FILE NAME FROM USER                        
!                                                                       
!    CALL LINE:    CALL GETFNM(FNAME,FNML,PROMPT,PROMTL)                
!                                                                       
!    INPUTS:       PROMPT = STRING TO POMPT USER WITH (C*60)            
!                  PROMTL = LENGTH OF PROMPT STRING (I)                 
!                                                                       
!    OUTPUTS:      FNAME = FILE NAME (C*60)                             
!                  FNML = FILE NAME LENGTH (I)                          
!***********************************************************************
!                                                                       
      SUBROUTINE GETFNM(FNAME,FNML,PROMPT,PROMTL) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER FNML, PROMTL 
      CHARACTER*60 FNAME, PROMPT 
!                                                                       
!                       Local variables                                 
!                                                                       
      INTEGER STDOUT, STDIN, I 
      DATA STDOUT/6/, STDIN/5/ 
!                                                                       
!                       Prompt user for file name                       
!                                                                       
      WRITE(STDOUT,'(T2,A,$)') PROMPT(1:PROMTL) 
      READ(STDIN,'(A)') FNAME 
!                                                                       
!                        Figure out input file name length              
      DO 10 I=1,20,1 
        IF(FNAME(I:I).EQ.' ') GOTO 20 
   10 END DO 
!                                                                       
   20 CONTINUE 
      FNML = I-1 
!                                                                       
!                                                                       
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         DMATRIX.FOR                                          
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATINV                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Inverts a N by N Matrix                              
!                                                                       
!    CALL LINE:    CALL MATINV(A,AINV,N)                                
!                                                                       
!    INPUTS:       A = Array to be inverted  (D)                        
!                  N = dimesion of arrays (I)                           
!                                                                       
!    OUTPUTS:      AINV = Inverse of A (D)                              
!                                                                       
!    CALLS :       Numerical recipes routines LUDCMP, LUBKSB            
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATINV(A,AINV,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N,N), AINV(N,N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER NPHYS, I, J 
      PARAMETER(NPHYS=10) 
      DOUBLE PRECISION TEMP(NPHYS,NPHYS), Y(NPHYS,NPHYS), D 
      INTEGER INDEX(NPHYS) 
!                                                                       
!                 Make a copy of the array, and initialize              
!                 the indentity matrix                                  
      DO 20 J=1,N,1 
        DO 10 I=1,N,1 
          Y(I,J) = 0.0 
          TEMP(I,J) = A(I,J) 
   10   CONTINUE 
        Y(J,J) = 1.0 
   20 END DO 
!                                                                       
!                                                                       
!                 LU decompose the matrix                               
      CALL LUDCMP(TEMP,N,NPHYS,INDEX,D) 
!                                                                       
!                                                                       
!                 Back substitute to get inverse                        
      DO 30 J=1,N,1 
        CALL LUBKSB(TEMP,N,NPHYS,INDEX,Y(1,J)) 
   30 END DO 
!                                                                       
!                                                                       
!                 Copy temporary array into the inverse array           
      DO 50 J=1,N,1 
        DO 40 I=1,N,1 
          AINV(I,J) = Y(I,J) 
   40   CONTINUE 
   50 END DO 
!                                                                       
!                                                                       
  999 RETURN 
      END                                           
!***********************************************************************
!                                                                       
!    MODULE:       MATADD                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Adds two N by M Matrices                             
!                                                                       
!    CALL LINE:    CALL MATINV(A,B,C,N,M)                               
!                                                                       
!    INPUTS:       A,B= Arrays to be added  (D)                         
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      C = Array containing A+B (D)                         
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATADD(A,B,C,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N,M 
      DOUBLE PRECISION A(N,M), B(N,M), C(N,M) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
          C(I,J) = A(I,J)+B(I,J) 
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATMUL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiplies two Matrices (LxM)x(MxN)                  
!                                                                       
!    CALL LINE:    CALL MATMUL(A,B,C,L,M,N)                             
!                                                                       
!    INPUTS:       A,B= Arrays to be added  (D)                         
!                  L,M,N = Dimensions of arrays (I)                     
!                                                                       
!    OUTPUTS:      C = Array containing A x B (D)                       
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATMUL(A,B,C,L,M,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER L, M, N 
      DOUBLE PRECISION A(L,M), B(M,N), C(L,N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J, K 
      DOUBLE PRECISION SUM 
!                                                                       
!                 Loop over all elements of the array                   
      DO 30 I=1,L,1 
        DO 20 J=1,N,1 
!                                                                       
!                 Initialize SUM for a new element                      
          SUM = 0.0 
!                 Calculate (i,j)th element                             
          DO 10 K=1,M,1 
            SUM = SUM+A(I,K)*B(K,J) 
   10     CONTINUE 
          C(I,J) = SUM 
!                                                                       
   20   CONTINUE 
   30 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!***********************************************************************
!                                                                       
!    MODULE:       MV_MUL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiplies a Matrix times a vector (NxN)x(N)         
!                                                                       
!    CALL LINE:    CALL MV_MUL(A,V,RV,N)                                
!                                                                       
!    INPUTS:       A = Array to be multiplied  (D)                      
!                  V = Vector to be multiplied (D)                      
!                  N = Dimensions of arrays & vector (I)                
!                                                                       
!    OUTPUTS:      RV = resultant vector (D)                            
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MV_MUL(A,V,RV,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N,N), V(N), RV(N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
      DOUBLE PRECISION SUM 
!                                                                       
!                 Loop over all elements of the array                   
      DO 20 I=1,N,1 
!                                                                       
!                 Initialize SUM for a new element                      
        SUM = 0.0 
!                 Calculate (i)th element                               
        DO 10 J=1,N,1 
          SUM = SUM+A(I,J)*V(J) 
   10   CONTINUE 
        RV(I) = SUM 
!                                                                       
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATSCL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiply a N by M Matrix by a scalar                 
!                                                                       
!    CALL LINE:    CALL MATSCL(A,SCALAR,B,N,M)                          
!                                                                       
!    INPUTS:       A = Array to be scaled  (D)                          
!                  SCALAR = Constant to multiply matrix by (D)          
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      B = Array containing SCALAR x A (D)                  
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATSCL(A,SCALAR,B,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N, M 
      DOUBLE PRECISION A(N,M), B(N,M), SCALAR 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
!                 Loop over all elements of the array                   
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
!                                                                       
          B(I,J) = SCALAR*A(I,J) 
!                                                                       
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATCOP                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Copy one N by M Matrix into another                  
!                                                                       
!    CALL LINE:    CALL MATCOP(A,B,N,M)                                 
!                                                                       
!    INPUTS:       A = Array to be copied  (D)                          
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      C = Array to be copied into (D)                      
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATCOP(A,B,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N,M 
      DOUBLE PRECISION A(N,M), B(N,M) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
          B(I,J) = A(I,J) 
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       VECCOP                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/10/90                                              
!                                                                       
!    PURPOSE:      Copy one vector on length N into another             
!                                                                       
!    CALL LINE:    CALL VECCOP(A,B,N)                                   
!                                                                       
!    INPUTS:       A = Array to be copied  (D)                          
!                  N = Number of rows in arrays (I)                     
!                                                                       
!    OUTPUTS:      B = Array to be copied into (D)                      
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE VECCOP(A,B,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N), B(N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I 
!                                                                       
      DO 10 I=1,N,1 
        B(I) = A(I) 
   10 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!                                                                       
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       LUDCMP                                               
!    TYPE:         SUBROUTINE                                           
!                                                                       
!***********************************************************************
!                                                                       
      SUBROUTINE LUDCMP(A,N,NP,INDX,D) 
!                                                                       
!               Make the default double precision                       
      implicit real*8(a-h,o-z) 
!                                                                       
!                                                                       
      PARAMETER (NMAX=100,TINY=1.0E-20) 
      DIMENSION A(NP,NP),INDX(N),VV(NMAX) 
      D=1. 
      DO 12 I=1,N 
        AAMAX=0. 
        DO 11 J=1,N 
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J)) 
   11   CONTINUE 
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.' 
        VV(I)=1./AAMAX 
   12 END DO 
      DO 19 J=1,N 
        IF (J.GT.1) THEN 
          DO 14 I=1,J-1 
            SUM=A(I,J) 
            IF (I.GT.1)THEN 
              DO 13 K=1,I-1 
                SUM=SUM-A(I,K)*A(K,J) 
   13         CONTINUE 
              A(I,J)=SUM 
            ENDIF 
   14     CONTINUE 
        ENDIF 
        AAMAX=0. 
        DO 16 I=J,N 
          SUM=A(I,J) 
          IF (J.GT.1)THEN 
            DO 15 K=1,J-1 
              SUM=SUM-A(I,K)*A(K,J) 
   15       CONTINUE 
            A(I,J)=SUM 
          ENDIF 
          DUM=VV(I)*ABS(SUM) 
          IF (DUM.GE.AAMAX) THEN 
            IMAX=I 
            AAMAX=DUM 
          ENDIF 
   16   CONTINUE 
        IF (J.NE.IMAX)THEN 
          DO 17 K=1,N 
            DUM=A(IMAX,K) 
            A(IMAX,K)=A(J,K) 
            A(J,K)=DUM 
   17     CONTINUE 
          D=-D 
          VV(IMAX)=VV(J) 
        ENDIF 
        INDX(J)=IMAX 
        IF(J.NE.N)THEN 
          IF(A(J,J).EQ.0.)A(J,J)=TINY 
          DUM=1./A(J,J) 
          DO 18 I=J+1,N 
            A(I,J)=A(I,J)*DUM 
   18     CONTINUE 
        ENDIF 
   19 END DO 
      IF(A(N,N).EQ.0.)A(N,N)=TINY 
      RETURN 
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       LUBKSB                                               
!    TYPE:         SUBROUTINE                                           
!                                                                       
!***********************************************************************
!                                                                       
!                                                                       
      SUBROUTINE LUBKSB(A,N,NP,INDX,B) 
!                                                                       
!               Make the default double precision                       
      implicit real*8(a-h,o-z) 
!                                                                       
!                                                                       
      DIMENSION A(NP,NP),INDX(N),B(N) 
      II=0 
      DO 12 I=1,N 
        LL=INDX(I) 
        SUM=B(LL) 
        B(LL)=B(I) 
        IF (II.NE.0)THEN 
          DO 11 J=II,I-1 
            SUM=SUM-A(I,J)*B(J) 
   11     CONTINUE 
        ELSE IF (SUM.NE.0.) THEN 
          II=I 
        ENDIF 
        B(I)=SUM 
   12 END DO 
      DO 14 I=N,1,-1 
        SUM=B(I) 
        IF(I.LT.N)THEN 
          DO 13 J=I+1,N 
            SUM=SUM-A(I,J)*B(J) 
   13     CONTINUE 
        ENDIF 
        B(I)=SUM/A(I,I) 
   14 END DO 
      RETURN 
      END                                           
                                                                        
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EL_EOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EL_EOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         2/12/91                                              
!                                                                       
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@SBAST1.SUNYSB.EDU or               
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!    PURPOSE:      The elctron and photon equation of state             
!                                                                       
!                                                                       
!    CALL LINE:    CALL EL_EOS(T,YE,BRYDNS)                             
!                                                                       
!    INPUTS:       T = TEMPERATURE                                      
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:      NONE                                                 
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EL_EOS.INC                                         
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EL_EOS(T,YE,BRYDNS) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION T, YE, BRYDNS 
!                                                                       
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                                                                       
!                           Plancks constant & speed of light           
      DOUBLE PRECISION HBAR, C 
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23) 
!                                                                       
!                           Pi and 1/3                                  
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23 
      PARAMETER(PI=3.1415927, PI2=9.8696044) 
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667) 
!                                                                       
!*****added by MBD                                                      
      double precision me, qcoef 
      me = 5.11d-1 
!                                                                       
!                    Leptons                                            
!                                                                       
!                    Electron number density                            
      NSUBE = BRYDNS*YE 
!                                                                       
!                    Coefficants for chemical potential                 
!                    and thermodynamics quantities                      
      QSUBE = 1.0/( 3.0*(PI**2)*((HBAR*C)**3) ) 
!                                                                       
      ACOEF = 0.5*NSUBE/QSUBE 
!                                                                       
!*****added by MBD new qcoef to include me term                         
      qcoef=((PI*T)**2)/3. - me*me/2. 
      BCOEF= (ACOEF**2+qcoef**3)**0.5 
!      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0)**0.5                      
!                                                                       
!      DBDT = (PI**6)*(T**5)/(9.0*BCOEF)                                
      DBDT = PI2*T*((PI2*T*T/3. -me*me/2.)**2)/BCOEF 
!                                                                       
      CCOEF = (ACOEF+BCOEF)**OVR3 
!                                                                       
!                                                                       
!                    Electron chemical potential                        
!      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF                            
      MUSUBE = CCOEF-qcoef/CCOEF 
!                                                                       
!*****added by MBD                                                      
!                                                                       
!      write(*,101) qcoef, ACOEF, BCOEF, CCOEF, MUSUBE                  
!                                                                       
!*****now do it again with me=0                                         
!                                                                       
!      qcoef=((PI*T)**2)/3.                                             
!      BCOEF= (ACOEF**2+qcoef**3)**0.5                                  
!      CCOEF = (ACOEF+BCOEF)**OVR3                                      
!      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF                            
!      write(*,101) qcoef, ACOEF, BCOEF, CCOEF, MUSUBE                  
!                                                                       
!101   format(1x,5(1pd11.4))                                            
!                                                                       
!                                                                       
!                                                                       
!                    Electron pressure for rel. case                    
!*****again me term added                                               
      EPRESS = 0.25*QSUBE*(MUSUBE**4+2.0*(PI*T*MUSUBE)**2+              &
     & 7.0*((PI*T)**4)/15.0 -3.0*(MUSUBE*me)**2 - 0.5*(PI*T*me)**2)     
!                                                                       
!                                                                       
!                    Electron internal energy per baryon                
!*****again me term added                                               
      EU = 0.75*QSUBE*(MUSUBE**4+2.0*(PI*MUSUBE*T)**2+                  &
     & 7.0*((PI*T)**4)/15.0 -(MUSUBE*me)**2 - 0.5*(PI*T*me)**2)         &
     & /BRYDNS                                                          
!                                                                       
!                                                                       
!                    Electron free energy per baryon                    
      FSUBE = ((MUSUBE*NSUBE)-EPRESS)/BRYDNS 
!                                                                       
!                    Electron entropy per baryon                        
      ES = QSUBE*(((PI*MUSUBE)**2)*T+7.0*(PI**4)*(T**3)/                &
     & 15.0 - 0.5*T*me*me)/BRYDNS                                       
!                                                                       
!                    Photons                                            
!                                                                       
!                    Photon pressure                                    
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3)) 
!                    Photon entropy per baryon                          
      PS = 4.0*PPRESS/(T*BRYDNS) 
!                                                                       
!                    Photon internal energy per baryon                  
      PU = 3.0*PPRESS/BRYDNS 
!                                                                       
!                    Photon free energy per baryon                      
      PF = PU-T*PS 
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                    Derivatives of chem. potential w.r.t. T,           
!                    BRYDNS, YE                                         
!                                                                       
!      DEMUDT = DBDT/(3.0*CCOEF**2)-OVR23*(PI**2)*T/CCOEF+              
!     1         DBDT*((PI*T)**2)/(9.0*CCOEF**4)                         
!                                                                       
!      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2)          
!                                                                       
      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2            &
     &         - me*me/2.)                                              
!                                                                       
      DEMUDT = -2.*DEMUDN*MUSUBE*T/(3.*YE*(HBAR*C)**3) 
!                                                                       
      DEMUDY = BRYDNS*DEMUDN/YE 
!                                                                       
!                                                                       
!                    Derivatives of pressure w.r.t. BRYDNS,YE,T         
!                                                                       
      DEPDN = BRYDNS*YE*DEMUDN 
!                                                                       
      DEPDY = BRYDNS*DEPDN/YE 
!                                                                       
      DEPDT = BRYDNS*(ES+YE*DEMUDT) 
!                                                                       
!                                                                       
!                    Derivatives of entropy w.r.t. T,BRYDNS,YE          
!                                                                       
      DESDT = ES/T+OVR23*(7.0*PI2*(T**2)/15.0+MUSUBE*T*DEMUDT)/         &
     &        (BRYDNS*(HBAR*C)**3)                                      
!                                                                       
      DESDN = -1.0*DEPDT/(BRYDNS**2) 
!                                                                       
      DESDY = 2.0*T*QSUBE*PI2*MUSUBE*DEMUDY/BRYDNS 
!                                                                       
!                                                                       
!                    Derivatives of internal energy w.r.t.              
!                    T,BRYDNS,YE                                        
      DEUDT = T*DESDT 
!                                                                       
      DEUDN = (YE*(MUSUBE-T*DEMUDT)-EU)/BRYDNS 
!                                                                       
      DEUDY = 3.0*QSUBE*((MUSUBE**3)+PI2*(T**2)*MUSUBE)*                &
     &        DEMUDY/BRYDNS                                             
!                                                                       
!                                                                       
!                               Photons                                 
!                                                                       
!                    Derivatives of photon pressure                     
      DPPDN = 0.0 
      DPPDT = BRYDNS*PS 
      DPPDY = 0.0 
!                                                                       
!                    Derivatives of photon entropy                      
      DPSDN = -PS/BRYDNS 
      DPSDT = 3.0*PS/T 
      DPSDY = 0.0 
!                                                                       
!                    Derivatives of internal energy                     
      DPUDN = -0.75*T*PS/BRYDNS 
      DPUDT = 3.0*PS 
      DPUDY = 0.0 
!                                                                       
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
                                                                        
!      Program to compute spline fits to fermi integrals                
!c  Must provide data file 14                                           
      subroutine initferm 
      implicit real*8(a-h,o-z) 
      parameter (n=201) 
      dimension f32(n),f12(n),fm12(n),eta(n),fr(n) 
      dimension f32a(n),f12a(n),fra(n),fia(n) 
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
!                                                                       
      open(14,file='fermi.atb',status='old') 
!                                                                       
      do 10 i=1,n 
       read(14,*)eta(i),f32(i),f12(i),fm12(i) 
   10  fr(i)=f12(i)/fm12(i) 
!                                                                       
      close(14,status='keep') 
!                                                                       
       call spline(eta,f12,n,f12a) 
!       write(*,1)f12a                                                  
       call spline(eta,f32,n,f32a) 
!       write(*,1)f32a                                                  
       call spline(eta,fr,n,fra) 
!       write(*,1)fra                                                   
       call spline(f12,eta,n,fia) 
!       write(*,1)fia                                                   
    1  format(1p8e10.3) 
       return 
      END                                           
                                                                        
      SUBROUTINE SPLINE(X,Y,N,Y2) 
!  Computes spline coefficients; Y(X) is input data; Y2 is output.      
        implicit real*8(a-h,o-z) 
        DIMENSION X(N),Y(N),Y2(N),U(500) 
        Y2(1)=0. 
        U(1)=0. 
        DO 11 I=2,N-1 
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1)) 
        P=SIG*Y2(I-1)+2. 
        Y2(I)=(SIG-1.)/P 
   11   U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))  &
         /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P 
        Y2(N)=0. 
        DO 12 K=N-1,1,-1 
   12   Y2(K)=Y2(K)*Y2(K+1)+U(K) 
        RETURN 
      END    

        SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,KLO,KHI) 
!     Computes spline fit of Y(X); YA(XA) is input data, Y2A are spline 
!  coefficents, klo and khi are running indices which bracket X.        
        implicit real*8(a-h,o-z) 
        DIMENSION XA(N),YA(N),Y2A(N) 
!c  Determine the bracketing indices                                    
        IF(KHI-KLO.GT.1)GOTO 1 
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
        KHI=KHI-1 
        KLO=KLO-1 
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
        KHI=KHI+2 
        KLO=KLO+2 
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
        KLO=1 
        KHI=N 
    1      IF(KHI-KLO.EQ.1) GOTO 2 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
        KHI=K 
        ELSE 
        KLO=K 
        ENDIF 
        GOTO 1 
    2      H=XA(KHI)-XA(KLO) 
        IF(H.EQ.0.) PAUSE 'BAD XA INPUT. ' 
!c  Compute spline fit.                                                 
        A=(XA(KHI)-X)/H 
        B=(X-XA(KLO))/H 
        Y=A*YA(KLO)+B*YA(KHI)+                              &
         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H**2/6. 
!	write(*,5)klo,khi,x,xa(klo),xa(khi),ya(klo),ya(khi)                   
!     > ,y2a(klo),y2a(khi),y                                            
    5      format(2i3,1p8e9.2) 
        RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       F_1_2                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    F_1_2(Y)      (1/2th Fermi Integral)                 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       1/2th Fermi Integral (DOUBLE PRECISION)              
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION F_1_2(Y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n), &
     & fia(n)                                                           
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.17677669&
     &,6.41500299D-02,.4D0,1.32934039D0,.33333333333d0,1,n/             
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,f12,f12a,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=A(6)*y*SQRT(y)*th*(5+(A(1)+(3*A(2)+7*X2*A(3))*X2)*x2) 
      GO TO 100 
   20 F0=DEXP(y) 
      F1=A(7)*th*F0*(2-(4*A(4)-(6*A(5)-.25*F0)*F0)*F0) 
  100 F_1_2=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       F_3_2                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    F_3_2(Y)      (3/2th Fermi Integral)                 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       3/2th Fermi Integral (DOUBLE PRECISION)              
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION F_3_2(y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n), &
     & fia(n)                                                           
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0&
     &,6.41500299D-02,.4D0,1.32934039D0,1,n/                            
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,f32,f32a,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=A(6)*SQRT(y)*(1./X2+A(1)-(A(2)+X2*A(3))*X2) 
      GO TO 100 
   20 F0=DEXP(y) 
      F1=A(7)*F0*(1.-(A(4)-(A(5)-.03125*F0)*F0)*F0) 
  100 F_3_2=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       FINV12                                               
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    FINV12(Y)      (Inverse of the 1/2th Fermi Integral) 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       Inverse of Fermi Integral (DOUBLE PRECISION)         
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION FINV12(y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION AI(8),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),        &
     & fra(n),f12a(n)                                                   
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA AI,klo,khi/-.822467032D0,-1.21761363D0,-9.16138616D0,        &
     &.398942281D0,.0732748216D0,-1.310707D0,1.12837917D0,              &
     &8.2810645D-3,1,n/                                                 
      if(y .gt. 109.695) goto 10 
      if(y .lt. 4.0234e-5) goto 20 
      call splint(f12,eta,fia,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=(1.5*y)**(.666666667) 
      X4=1./(X2*X2) 
      F1=X2*(1.+(AI(1)+(AI(2)+AI(3)*X4)*X4)*X4) 
      GO TO 100 
   20 F1=LOG(AI(7)*MAX(y,1.D-20)*(1.+(AI(4)+(AI(5)+AI(8)*y)*y)*y)) 
  100 finv12=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       FHALF                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         3/10/90                                              
!                                                                       
!    CALL LINE:    FHALF(Y)                                             
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       Ratio of 1/2th Fermi Integral to the -1/2th Fermi    
!                  Integral (DOUBLE PRECISION)                          
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION FHALF(y) 
      IMPLICIT real*8(a-h,o-z) 
      parameter(n=201) 
      DIMENSION A(7),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),         &
     & fra(n),f12a(n)                                                   
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.17677669&
     &,6.41500299D-02,.4D0,1.32934039D0,.3333333333333d0,1,n/           
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,fr,fra,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=y*th*(1.+(.2*A(1)+(.6*A(2)+1.4*X2*A(3))*X2)*x2)                &
     & /(1.-(.2*th*a(1)-(a(2)-4.2*x2*a(3))*x2)*x2)                      
      GO TO 100 
   20 F0=EXP(y) 
      F1=(1.-(2*a(4)-(3*a(5)-.125*f0)*f0)*f0)/                          &
     & (2.-(8*a(4)-(18*a(5)-f0)*f0)*f0)                                 
  100 FHALF=F1 
  999 RETURN 
      END                                           
