      program lahyc 
!*******************************************************                
!                                                      *                
!  This is a Lagrangian 1D hydrodynamics code.         *                
!  adds particles                                      *                
!  this is the main driver of the integration.         *                
!                                                      *                
!*******************************************************                
!      
      implicit double precision (a-h,o-z) 
!                                                                       
!--ntstep counts the number of timesteps                                
      integer jtrape,jtrapb,jtrapx,ntstep
      character*128 mlmodel_name
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /carac/ deltam(idim), abar(idim) 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /cgas / gamma 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /ener2/ tkin, tterm 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /epcap/ betafac, c2cu, c3cu 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /freez/ ufreez(idim) 
      common /numb / ncell, ncell1 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)                  
      common /propt/ dtime,tmax 
      common /shock/ cq,cl 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common /tempe/ temp(idim) 
      common /therm/ xmu(idim) 
      common /timej / time, dt
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x

      logical trapnue, trapnueb, trapnux, print_nuloss 
      common /trap / trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /typef/ iextf, ieos 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
      common /uswest/ usltemp, uslrho, uslu, uslp, u2slu
      common /mlmod/ mlmodel_name
      common /idump/ idump
!                                                                       
!--read initial condition   
!
      call readini      
!                                                                       
!--define code's units                                                  
!                                                                       
      call unit 
!                                                                       
!--set all quantities needed for run                                    
!                                                                       
      call preset 
!                                                                       
!--main loop to save the dump file                                      
!                                                                                                          
!--ntstep is the number of timesteps                                          
      ntstep = 1 
          
      do while (time.lt.tmax) 
!42   continue
        !--idump is a counter for dump files         
        idump = idump + 1              
!                                                                    
!--step runs the runge-kutta operator                                   
!                                                                       
        call step(ntstep) 
!                                                                       
!--produce output if desired                                            
!                                                                       
        lu=61 
        call printout(lu) 
        print*,' '    
        print*,'****************' 
        print*,' idump        = ', idump
        print*,' savetime (s) = ', time*10
        print*,'****************' 
        print*,' ' 
!      if (time.gt.4.0) dtime = .0001                                   
!                                                                       
!--check if simulation is finished                                      
!                                                                       
        !print *,time,tmax                                              
        !if(time.gt.tmax) go to 90                                      
!                                                                       
        !go to 42                                                              
      end do 
!                                                                       
      !90 call printout(lu)                                             
      !call printout(lu) 
      stop 
      END                                           
!
! ============================ Subroutines ============================
!      
      subroutine hydro(time,ncell,x,v,                                  &
     &           u,rho,ye,f,du,dye,q,                                   &
     !fmix,                              &
     &           ynue,ynueb,ynux,dynue,dynueb,dynux,                    &
     &           unue,unueb,unux,dunue,dunueb,dunux,                    &
     &           print_nuloss,ntstep)                                          
!                                                                       
!****************************************************************       
!                                                               *       
!  this subroutine advances the system of hydro equations by    *       
!  one time step.                                               *       
!                                                               *       
!  Description of the variables :                               *       
!  ------------------------------                               *       
!                                                               *       
!  ncell        :  number of cells                              *       
!  ncell1       :  number of cell edges                         *       
!  rho          :  density (cell centered)                      *       
!  pr           :  pressure (cell centered)                     *       
!  u            :  specific internal energy (cell centered)     *       
!  q            :  artificial viscous stress (cell centered)    *       
!  deltam       :  mass of the cell (cell centered)             *       
!  prold,qold,rhold  : pressure artificial viscous stress and   *       
!                     density at previous time step             *       
!  x            :  cell boundaries (edge centered)              *       
!  v            :  velocity (edge centered)                     *       
!  vold         :  velocity at previous time step               *       
!  gamma        :  ratio of specific heat                       *       
!  cq, cl       :  quadratic and linear artificial viscous      *       
!                  stress coefficients                          *       
!  tkin, uint   :  kinetic and specific internal energy         *       
!  time         :  current time                                 *       
!  dt, dtold    :  current and previous time step               *       
!  tmax         :  maximum time for the simulation              *       
!                                                               *       
!****************************************************************       
!                                                                       
!  cell                                                                 
! center          1         2         3                                 
!------------|---------|---------|---------|---                         
!  cell      0         1         2         3                            
!  edge                                                                 
!                    
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx, mlin_grid_size
      character*128 mlmodel_name
      logical post_bounce

      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (tiny=-1e-5) 
!                                                                       
      dimension x(0:idim), v(0:idim), f(0:idim) 
      dimension fmix(idim) 
      dimension u(idim), rho(idim), ye(idim) 
      dimension dye(idim), du(idim), q(idim) 
      dimension ynue(idim),ynueb(idim),ynux(idim) 
      dimension unue(idim),unueb(idim),unux(idim) 
      dimension dynue(idim),dynueb(idim),dynux(idim) 
      dimension dunue(idim),dunueb(idim),dunux(idim) 
!                                                                       
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /carac/ deltam(idim), abar(idim) 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /ener1/ dq(idim), dunu(idim) 
      common /ener2/ tkin, tterm 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /epcap/ betafac, c2cu, c3cu 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /freez/ ufreez(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /shock/ cq,cl 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common /tempe/ temp(idim) 
      common /therm/ xmu(idim) 
      logical trapnue, trapnueb, trapnux, print_nuloss 
      common /trap / trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /typef/ iextf, ieos 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
      common /mlmod/ mlmodel_name
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
      common /interp/ mlin_grid_size
      common /mlout/ pr_turb(idim1)
!      common /nuout/ rlumnue, rlumnueb, rlumnux,                       
!     1               enue, enueb, enux, e2nue, e2nueb, e2nux           
!                                        
!--compute density                                                      
! ---------------------------------------------------------             
!    
      call density(ncell,x,rho) 
!                                                                      
!                                                                       
!--compute thermodynamical properties                                   
!  ----------------------------------                                   
!  (this allows the use of simple eos)                                  
!  (ieos=3 or 4 calls for the sophisticated eos's)                      
!                               
      !do i=30,40
      !    print*, i, ye(i)
      !enddo
      !print *, '[time] ',time,ieos                                     
      if(ieos.eq.1)call eospg(ncell,rho,u) 
      if(ieos.eq.2)call eospgr(ncell,rho,u) 
      if(ieos.eq.3.or.ieos.eq.4)call eos3(ncell,rho,u,ye) 
      if(ieos.eq.5)call eos5(ncell,rho,u,ye) 
!                                                                       
!--do gravity                                                           
!  --------------------------------------------------------             
!                                                                       
      call gravity(ncell,deltam,x,f) 
!                                                                       
!                                                                       
!--neutrino physics                                                     
!  ----------------                                                     
      yemn=1e20 
      yemx=-1e20 
      ynuemx=-1e20 
      ynuebmx=-1e20 
      ynuxmx=-1e20 
      do k=1,ncell 
         yemn=dmin1(ye(k),yemn) 
         yemx=dmax1(ye(k),yemx) 
         ynuemx=dmax1(ynue(k),ynuemx) 
         ynuebmx=dmax1(ynueb(k),ynuebmx) 
         ynuxmx=dmax1(ynux(k),ynuxmx) 
         if (ye(k).lt.0.0) print *,'ye<0, k=',k,ye(k) 
!         if (ye(k).gt.0.51) print *,'ye>0.51, k=',k,ye(k)              
         if (ynue(k).lt.tiny) print *,'ynue<0, k=',k,ynue(k),            &
     &                               trapnue(k),ebetaeq(k)              
                                                                        
         if (ynux(k).lt.tiny) print*,'ynux<0, k=',k,ynux(k) 
      enddo 
!      write(*,200)'hydro: yemn,ynuemx,ynuebmx,ynuxmx',                 
!     1                   yemn,ynuemx,ynuebmx,ynuxmx                    
  200 format(A33,4(1x,1pe10.3)) 
!                                                                       
!-- initialize                                                          
!--skip neutrino physics                                                
!      goto 90                                                          
!                                                                       
      call nuinit(ncell,rho,x,ye,dye,                                   &
     &            ynue,ynueb,ynux,dynue,dynueb,dynux,                   &
     &            unue,unueb,unux,dunue,dunueb,dunux)                   
!                                                                                                           
!-- nu contribution to pressure                                         
!                                                                       
      call nupress(ncell,rho,unue,unueb,unux) 
!                                                                       
!--e+/e- capture                                                        
!                                                                       
      call nuecap(ncell,rho,ye,dye,dynue,dynueb,dunue,dunueb) 
!                                                                       
!--plasma/pair neutrino emission processes                              
!                                                                       
      call nupp(ncell,rho,ye,dynue,dynueb,dynux,                        &
     &          dunue,dunueb,dunux)                                     
!                                                                       
!--neutrino/anti neutrino conversion                                    
!                                                                       
      call nuconv(ncell,x,rho,ynue,ynueb,ynux,                          &
     &     dynue,dynueb,dynux,dunue,dunueb,dunux)                       
!                                                                       
!--neutrino/anti-neutrino annihilation                                  
!                                                                       
!      call nuann(ncell,x,rho,ynue,ynueb,ynux,                          
!     $          dynue,dynueb,dynux,dunue,dunueb,dunux)                 
!                                                                       
!--neutrino depletion from thin regions                                 
!                                                                       
      call nusphere(ncell,x,                                            &
     &            ynue,ynueb,ynux,dynue,dynueb,dynux,                   &
     &            unue,unueb,unux,dunue,dunueb,dunux)                   
!                                                                       
!                                                                       
!--fold-in neutrino background luminosity from the core                 
!--and normalize energy sums                                            
!                                                                       
      call nulum(print_nuloss) 
!                                                                       
!--neutrino absorption                                                  
!                                                                       
      call nuabs(ncell,rho,x,dye,ynue,ynueb,                            &
     &          dynue,dynueb,dunue,dunueb)                              
!                                                                       
!--neutrino/electron scattering                                         
!                                                                       
      call nuscat(ncell,rho,x,ynue,ynueb,ynux,                          &
     &            dunue,dunueb,dunux)                                   
!                                                                       
!--do the beta equilibrium cases                                        
!                                                                       
      call nubeta(ncell,x,rho,ye,dye,ynue,ynueb,unue,unueb,             &
     &            dynue,dynueb,dunue,dunueb)                            
!                                                                       
   90 continue 
!                                                                       
!--compute turbulence parameters                                        
!                                                                       
      !call turbulence(ncell,x,f,q,v,rho,fmix) 
!                                                                       
!--compute q values                                                     
!                                                                       
      call artvis(ncell,x,rho,v,q)                                                                     
!     
      !post_bounce = .true.
      if (post_bounce.eqv..true.) then
        !--calculate PNS & shock radii, only in post-bounce stage
        call shock_radius(ncell,x,v,print_nuloss)
        call pns_radius(ncell,x,rho,print_nuloss)
        
        !--turbulence contribution to pressure via ML in post-bounce regime
        if (mlmodel_name == 'None') then
            pr_turb = 0
        else
            call turbpress(ncell,rho,x,v,temp)
        endif
      else
          !--check if bounced
          call bounce(ntstep)
      endif
!
!--compute forces on the particles                                      
      call forces(ncell,x,f,q,v,rho) 
!                                                                       
!--flux limited diffusion                                               
!--skip neutrino                                                        
      call nudiff(ncell,x,rho,time,ye,                                  &
     &    ynue,ynueb,ynux,dynue,dynueb,dynux,                           &
     &    dunue,dunueb,dunux)                                           
!                                                                       
!--compute energy derivative                                            
!                                                          
      call energ(ncell,x,v,dye,du,rho) 
!                                                                       
!--compute neutrino pressure work and change of <E>s                    
!                                                                       
!--skip neutrino                                                        
      call nuwork(ncell,x,v,rho,                                        &
     &            unue,unueb,unux,dunue,dunueb,dunux)                   
!                                                                       
   99 return 
      END                                           
!
!     
      subroutine turbpress(ncell,rho,x,v,temp) 
!*********************************************************               
!                                                        *               
! This subroutine passes the grid to a trained           *
! PyTorch model that predicts pressure due to turbulence *                             
!                                                        *               
!*********************************************************                
!
      use pytorch
      use data_functions
      implicit double precision (a-h,o-z) 
!          
      character*128 :: mlmodel_name
      integer mlin_grid_size
      real scale_pr, scale_pr_relative
      double precision pr_turb
!
      parameter(idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension rho(idim), temp(idim) 
      dimension x(0:idim), v(0:idim)
      dimension unue(idim),unueb(idim),unux(idim)       
!                                                                       
      logical trapnue, trapnueb, trapnux 
      
      common /trap/ trapnue(idim),trapnueb(idim),trapnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /eosnu/ prnu(idim1) 
      common /mlmod/ mlmodel_name               
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x       
      common /interp/ mlin_grid_size
      common /mlout/ pr_turb(idim1)
                                    
      ! The tensor shape is exactly backwards from python: (Length,Channels,N batches)
      real(real32) :: input(mlin_grid_size, 4, 1)
      real(real32), allocatable :: output_h(:,:,:)    
      double precision pr_relative(idim1)
      double precision interp_x(mlin_grid_size)         
      
      ! Scale Pressure to fit into single precission (taken from ML training)
      scale_pr = 1.d-8
      scale_pr_relative = 1.d-8

      input(:,1,1) = interpolate(x(1:),v(1:),ncell,int(pns_ind),int(shock_ind),mlin_grid_size)
      input(:,2,1) = interpolate(x(1:),rho,ncell,int(pns_ind),int(shock_ind),mlin_grid_size)
      input(:,3,1) = interpolate(x(1:),pr(1:),ncell,int(pns_ind),int(shock_ind),mlin_grid_size)*scale_pr
      input(:,4,1) = interpolate(x(1:),temp,ncell,int(pns_ind),int(shock_ind),mlin_grid_size)

      ! ML model (reverse the order for fortran data)
      ! Input: ['u1','rho','Pgas', 'T'] (1, 4, 200)
      ! Output: [pr_relative = P_turb/P_gas] (1, 1, 200)

      output_h = mlmodel(input, trim(mlmodel_name))      
      
      ! Re-shape output into code-shape mlout:
      pr_turb = 0   
      pr_relative = 0
      interp_x = linspace(x(int(pns_ind)), x(int(shock_ind)), mlin_grid_size)   
      pr_relative(pns_ind:shock_ind) = interpolate(DBLE(interp_x),DBLE(output_h(:,1,1)),      &
                                                   mlin_grid_size,1,mlin_grid_size,           &
                                                   int(shock_ind-pns_ind))*scale_pr_relative
      pr_turb = pr_relative*pr
      print*, 'ML prediction'      
      print*, size(output_h), shape(output_h), shock_ind-pns_ind
      print*, '-- pr_relative --'
      print*, pr_relative(shock_ind-5:shock_ind)
      print*, '-- pr_turb --'
      print*, pr_turb(shock_ind-5:shock_ind)
      call exit(0)
      return 
      END       
!
!
      subroutine shock_radius(ncell,x,v,print_nuloss)
!***********************************************************            
!                                                          *            
!  This subroutine identifies the shock radius             *            
!                                                          *            
!***********************************************************            
!            
      implicit double precision (a-h,o-z) 

      common /rshock/ shock_ind, shock_x    
      dimension x(0:ncell),v(0:ncell) 
      real :: initial_v, old_max_v
      real, dimension(ncell) :: v_old
      integer :: i,j, ind
      logical print_nuloss
!                              
      !initial_v = v(size(v))      

      !do i=size(v), 1, -1
      !    if (v(i).le.(initial_v-(initial_v-minval(v))*0.1)) then
      !        old_max_v = v(i)
      !        do j=i-1,0,-1
      !            if (v(j) .le. old_max_v) then
      !                old_max_v = v(j)
      !            else
      !                shock_x = x(j)
      !                shock_ind = j
      !                EXIT
      !            end if 
      !        end do
      !        EXIT
      !    end if
      !end do
      
      shock_ind = minloc(v, dim=1)
      shock_x = x(shock_ind)

511   format(A,1p,I5,A,E10.3) 
      if (print_nuloss .eqv. .true.) then      
          write(*,511)'[ shock radius (i, km) ]', int(shock_ind),               &
          '                    ', 1.d4*shock_x          
      end if 
      
      return
      end
!
      subroutine pns_radius(ncell,x,rho,print_nuloss)
!***********************************************************            
!                                                          *            
!  This subroutine identifies the Proto-Neutron Star radius*            
!                                                          *            
!***********************************************************            
!            
      implicit double precision (a-h,o-z) 

      common /pns/ pns_ind, pns_x   
      dimension x(0:ncell)
      dimension rho(ncell) 
      real :: rho_threshold      
      integer :: i
      logical print_nuloss
!      
!--g/cm^3 / unit conversion
      rho_threshold = 1.d13/2.d6
!
      do i=size(rho), 1, -1        
          if (rho(i) .ge. rho_threshold) then                        
              pns_x = x(i)
              pns_ind = i
              EXIT
          end if 
      end do      

511   format(A,1p,I5,A,E10.3) 
      if (print_nuloss .eqv. .true.) then      
          write(*,511)'[   PNS radius (i, km) ]', int(pns_ind),               &
          '                    ', 1.d4*pns_x          
      end if 
      
      return
      end
!
!
      subroutine bounce(ntstep)
!***********************************************************            
!                                                          *            
!  This subroutine identifies if the bounce has occured    *            
!                                                          *            
!***********************************************************            
!            
      implicit double precision (a-h,o-z) 

      logical post_bounce

      parameter (idim=10000) 

      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
           &               enue, enueb, enux, e2nue, e2nueb, e2nux
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
      common /timej / time, dt
      
      post_bounce = .false.
      
      if (ntstep==1) then
          rlumnue_max = rlumnue
      endif
      
      if (rlumnue>(2*rlumnue_max)) then
          post_bounce = .true.
          bounce_ntstep = ntstep
          bounce_time = time
      else if (rlumnue > rlumnue_max) then
          rlumnue_max = rlumnue
      endif                  
      
      end
!
      subroutine artvis(ncell,x,rho,v,q) 
!***********************************************************            
!                                                          *            
!  This subroutine updates the q-value for artificial      *            
!  viscosity.                                              *            
!                                                          *            
!***********************************************************            
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension x(0:idim),v(0:idim) 
      dimension rho(idim),q(idim) 
!                                                                       
      common /ener1/ dq(idim), dunu(idim) 
      common /shock/ cq,cl 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /carac/ deltam(idim), abar(idim) 
!                                                                       
      data pi4/12.56637d0/ 
!                                                                       
!--update q value                                                       
!                                                                       
      do kp05=1,ncell 
         k=kp05 - 1 
         k1=kp05 
         q(kp05)=0.d0 
         akp1=pi4*x(k1)*x(k1) 
         ak=pi4*x(k)*x(k) 
         akp05=0.5d0*(ak+akp1) 
         gradv=v(k1)*(3.d0*akp05-akp1) - v(k)*(3.d0*akp05-ak) 
         if(gradv.lt.0.d0)then 
!                                                                       
!--quadratic term                                                       
!                                                                       
            dv=v(k1) - v(k) 
            alpha=0.5d0*cq*cq*rho(kp05)*dabs(dv)/akp05 
            q(kp05)=-alpha*gradv 
!                                                                       
!--linear term                                                          
!                                                                       
            cs=vsound(kp05) 
            alphal=0.5d0*cl*rho(kp05)*cs/akp05 
            q(kp05)=q(kp05) - alphal*gradv 
         end if 
         dq(kp05)=-q(kp05)*gradv/deltam(kp05) 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine coulomb(rhok,zbar,ye,ucoul,pcoul) 
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
      parameter(ufac=-0.214d0) 
      parameter(pfac=-0.0714d0) 
!                                                                       
      rho13=rhok**0.333333333333d0 
      rho43=rho13*rhok 
      ye2=ye*ye 
      y43z23=(ye2*zbar)**0.66666666666d0 
      ucoul=ufac*rho13*y43z23 
      pcoul=pfac*rho43*y43z23 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine density(ncell,x,rho) 
!****************************************************************       
!                                                               *       
!  This subroutine calculates the density using the continuity  *       
!  equation.                                                    *       
!                                                               *       
!****************************************************************       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
!                                                                       
      dimension x(0:idim),rho(idim) 
      common /carac/ deltam(idim), abar(idim) 
!                                                                       
      data pi4/12.566371/ 
      integer :: i
!                                                                       
!--update density                                                       
!                                                                       
      do kp05=1,ncell 
         k1=kp05 
         k=kp05 - 1 
         rho(kp05)=3.d0*deltam(kp05)/                                   &
     &              (pi4*(x(k1)*x(k1)*x(k1)-x(k)*x(k)*x(k)))            
      enddo
      
      return 
      END                                           
!                                                                       
      subroutine energ(ncell,x,v,dye,du,rho) 
!************************************************************           
!                                                           *           
!  this routine computes the change in internal energy      *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (avokb=6.02e23*1.381e-16) 
!                                                                       
      dimension x(0:idim),v(0:idim) 
      dimension du(idim), dye(idim), rho(idim) 
!                                                                       
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /tempe/ temp(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /typef/ iextf, ieos 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /ener1/ dq(idim), dunu(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /carac/ deltam(idim), abar(idim) 
      !common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
      double precision dupp 
!                                                                       
      data pi4/12.566371/ 
!                                                                       
!--compute change in specific internal energy                           
!                                                                       
!--entropy conversion factor                                            
      sfac=avokb*utemp/uergg 
!                                                                       
      do kp05=1,ncell 
         k=kp05-1 
         k1=kp05 
         !vturbk=dsqrt(vturb2(kp05)) 
         akp1=pi4*x(k1)*x(k1) 
         akp=pi4*x(k)*x(k) 
         !pdv=(pr(kp05)+rho(kp05)*vturb2(kp05))*                         &
     !&        (akp1*v(k1)-akp*v(k))/deltam(kp05)                        
         pdv=pr(kp05)*(akp1*v(k1)-akp*v(k))/deltam(kp05)
         
!                                                                       
         dupp=0.0 
         if (temp(kp05).lt.6.and.rho(kp05).lt.1000.) then 
            dupp=9.96d5*temp(kp05)**9/rho(kp05) 
         end if 
         if (ifleos(kp05).ne.3) then 
!-- we subtract the energy added to the neutrino field                  
! no shock heating                                                      
            !du(kp05)=-pdv+0.5*dq(kp05)-dunu(kp05)-dupp+                 &
     !&           rho(kp05)*vturbk**3/dmix(kp05)                         
            du(kp05)=-pdv+0.5*dq(kp05)-dunu(kp05)-dupp     
            if (rho(kp05+1).gt.0.1.and.kp05.lt.1) then 
               xp05=.5d0*(x(k)+x(k1)) 
               xp15=.5d0*(x(k1)+x(k1+1)) 
               if (temp(kp05)*xp05.gt.1.2d0*temp(kp05+1)*xp15) then 
                  dubef=du(kp05) 
                  du(kp05)=min(1.2d0*du(kp05),-1.d-2*du(kp05)) 
                  print *, 'Artificial Cooling', dubef, du(kp05) 
               end if 
            end if 
!            du(kp05)=-pdv - dunu(kp05)                                 
         else 
            if (ieos.eq.3) then 
!--the energy is the variable of state:                                 
            du(kp05)=-pdv + 0.5*dq(kp05) - dunu(kp05) 
            if (rho(kp05+1).gt.0.1.and.kp05.lt.1) then 
               xp05=.5d0*(x(k)+x(k1)) 
               xp15=.5d0*(x(k1)+x(k1+1)) 
               if (temp(kp05)*xp05.gt.1.2d0*temp(kp05+1)*xp15) then 
                  dubef=du(kp05) 
                  du(kp05)=min(1.2d0*du(kp05),-1.d-2*du(kp05)) 
                  print *, 'Artificial Cooling', dubef, du(kp05) 
               end if 
            end if 
!            du(kp05)=-pdv  - dunu(kp05)                                
!                                                                       
            elseif (ieos.eq.4) then 
!--the entropy is the variable of state:                                
               du(kp05)=(0.5*dq(kp05) - dunu(kp05) +                    &
     &                sfac*dye(kp05)*(xmuhat(kp05)-                     &
     &                xmue(kp05)))/temp(kp05)                           
               if (rho(kp05+1).gt.0.1.and.kp05.lt.1) then 
                  xp05=.5d0*(x(k)+x(k1)) 
                  xp15=.5d0*(x(k1)+x(k1+1)) 
                  if (temp(kp05)*xp05.gt.1.2d0*temp(kp05+1)*xp15) then 
                     dubef=du(kp05) 
                     du(kp05)=min(1.2d0*du(kp05),-1.d-2*du(kp05)) 
                     print *, 'Artificial Cooling', dubef, du(kp05) 
                  end if 
               end if 
!     du(kp05)=(-dunu(kp05) +                                           
!     1                sfac*dye(kp05)*(xmuhat(kp05)-                    
!     2                xmue(kp05)))/temp(kp05)                          
!               if (kp05.lt.5) then                                     
!                  ent=sfac*dye(kp05)*(xmuhat(kp05)-                    
!     $                 xmue(kp05))                                     
!                  print *, kp05,du(kp05),dq(kp05),dunu(kp05),ent       
!               endif                                                   
            endif 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine eosflg(ncell,rho,ye,u,dye,du) 
!****************************************************************       
!                                                                       
! this subroutine determines what kind of eos to use:                   
!       eosflg = 1: freeze-out, just Ocean's eos + Coul corr.           
!       eosflg = 2: NSE with Raph's routines, + Ocean eos + Coul        
!       eosflg = 3: Swesty's eos   
!                                                                       
! and if e-neutrinos or x-neutrinos are trapped                         
!                                                                       
!**************************************************************         
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (iqn=17) 
      parameter (avokb=6.02e23*1.381e-16) 
!                                                                                     
      real ycc,yccave       
      double precision umass 
      double precision rhok, dens, tempk, yek, xpk, xnk, xak, xhk, yehk 
      double precision zbark, abark, abar2, ubind, dubind 
      double precision inpvar(4), xmuh 
      double precision brydns,pprev,psl,usl,dusl,gamsl,etak 
      double precision ucoul, pcoul 
      double precision pel,eel,sel,ptot,etot,stot,dpt,det,dpd,ded 
!                                                                       
      dimension rho(idim), ye(idim), u(idim) 
      dimension dye(idim), du(idim) 
!                                                                       
      common /ceos / amas(iqn), znum(iqn) 
      common /cc   / ycc(idim,iqn), yccave(iqn) 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /carac/ deltam(idim), abar(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /freez/ ufreez(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /tempe/ temp(idim) 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu 
      common /typef/ iextf, ieos 
!                                                                       
      data pi4/12.56637d0/ 
!                                                                       
      kount12=0 
      kount21=0 
      kount23=0 
      kount32=0 
      rhomax=0.0 
      t9freeze=t9nse-1.5 
!                                                                       
!--entropy conversion factor                                            
      sfac=avokb*utemp/uergg 
!                                                                       
      do k=1,ncell 
         tempk=temp(k) 
         rhok=rho(k)*udens 
         yek=ye(k) 
         rhomax=dmax1(rhomax,rhok) 
         if(ifleos(k).eq.1.and.temp(k).gt.t9nse) then 
            call nsestart(tempk,rhok,yek,xpk,xnk) 
            xp(k)=xpk 
            xn(k)=xnk 
            abar(k)=abark 
!                                                                       
!--for NSE, add in the nuclear component to the thermal                 
!--energy to get the total available internal energy                    
!           
            print *,'For NSE, nuclear component was added to E_thermal to get E_total,'
            print *,'change occured at cell',k,u(k),ufreez(k)           &
     &           ,ifleos(k)                                             
            u(k)=u(k)+ufreez(k) 
            ufreez(k)=0.0 
            ifleos(k)=2 
            kount12=kount12+1 
         endif 
         if(ifleos(k).eq.2.and.temp(k).le.t9freeze) then 
            ifleos(k)=1 
            kount21=kount21+1 
            xpk=xp(k) 
            xnk=xn(k) 
            call nsetemp(k,tempk,rhok,yek,tempk,xpk,xnk,                &
     &                   xak,xhk,yehk,zbark,abark,ubind,dubind)         
!                                                                       
!--for freeze-out, remove the nuclear component from the totalc--energy 
!                            
            print *,'For freeze-out, nuclear component was removed from E_total'
            print *,'change occured at cell',k,u(k),-ubind/uergg        &
     &           ,ifleos(k)                                             
            u(k)=u(k)-ubind/uergg 
            ufreez(k)=ubind/uergg 
            abar(k)=abark 
            dyccm=dmod(abark,4.0d0) 
            iycc=int(abark/4.d0) 
            if (iycc.le.2) then 
               dyccm=dyccm/8.d0 
               iycc=4 
            elseif (iycc.ge.15) then 
               dyccm=0 
               iycc=17 
               write (21,*) 'beyond the network', k, abark 
            else 
               dyccm=dyccm/4.d0 
               iycc=iycc+2 
            end if 
            do iq=1,17 
               ycc(k,iq)=0. 
            end do 
            ycc(k,1)=real(yek) 
            ycc(k,2)=real(xpk) 
            ycc(k,3)=real(xnk) 
            ycc(k,iycc)=(1.0-real(dyccm))/(amas(iycc))**2 
            ycc(k,iycc+1)=real(dyccm)/(amas(iycc))**2 
         elseif(ifleos(k).ne.3.and.rho(k).gt.rhoswe) then 
!ifleos(k).eq.2.and.rho(k).gt.rhoswe) then                              
            ifleos(k)=3 
            kount23=kount23+1 
!--make call to sleos to get the entropy or intenal energy              
!--at present rho, ye, T                                                
            brydns=rho(k)*uslrho 
            pprev=xpf(k) 
            inpvar(1)=tempk*usltemp 
            inpvar(2)=pvar2(k) 
            inpvar(3)=pvar3(k) 
            inpvar(4)=pvar4(k) 
!            if (k.gt.20.and.k.lt.30) then                              
!               print *, k, rho(k), tempk, u(k)                         
!               print *, k, inpvar(1), yek, brydns                      
!            end if                                                     
            call slwrap(k,inpvar,yek,brydns,pprev,                      &
     &           psl,usl,dusl,gamsl,etak,xpk,xnk,xak,xhk,               &
     &           yehk,abark,xmuh,stot)                                  
            xpf(k)=pprev 
            pvar2(k)=inpvar(2) 
            pvar3(k)=inpvar(3) 
            pvar4(k)=inpvar(4) 
            u(k)=usl/uslu 
            xp(k)=xpk 
            xn(k)=xnk 
            eta(k)=etak 
            xmuek=etak*temp(k) 
            xmuhk=xmuh/utmev 
!--recalculate du if necessary                                          
            if (ieos.eq.4) then 
               du(k)=(0.5*dq(k) - dunu(k) +                             &
     &                 sfac*dye(k)*(xmuhk-xmuek))/                      &
     &                 temp(k)   
            print *,'du was recalculated,'
            print *,'change occured at cell',k,u(k),ifleos(k) 
            endif 
         elseif(ifleos(k).eq.3.and.rho(k).lt.rhoswe) then 
!-- switch back to internal energy variable of state                    
            call nsestart(tempk,rhok,yek,xpk,xnk) 
            call nsetemp(k,tempk,rhok,yek,tempk,xpk,xnk,                &
     &                   xak,xhk,yehk,zbark,abark,ubind,dubind)         
            dens=rho(k)*uorho1 
            abar2=zbark/yek 
            call nados(tempk,dens,zbark,abar2,pel,eel,sel,              &
     &                  ptot,etot,stot,dpt,det,dpd,ded,gamsl,etak)      
            dens=rho(k) 
            call coulomb(dens,zbark,yek,ucoul,pcoul) 
            u(k)=ucoul+ubind/uergg+etot/uou1 
            xp(k)=xpk 
            xn(k)=xnk 
            eta(k)=etak 
            ifleos(k)=2 
            kount32=kount32+1 
!--recalculate du if needed                                             
            if (ieos.eq.4) then 
               km1=k-1 
               akp1=pi4*x(k)*x(k) 
               akp=pi4*x(km1)*x(km1) 
               pdv=pr(k)*(akp1*v(k)-akp*v(km1)) 
!-- we subtract the energy added to the neutrino field                  
               du(k)=-pdv/deltam(k) + 0.5*dq(k) - dunu(k) 
            endif 
            !print *,'change occured at cell',k,ifleos(k)
         endif 
!                                                                       
      enddo 
!      write(*,100) rhomax,kount12,kount21,kount23,kount32              
!  100 format('eosflg: rhomax, 1->2, 2->1, 2->3, 3->2',1pe10.2,4(1x,I4))
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine eospg(ncell,rho,u) 
!c************************************************************          
!                                                           *           
!  This subroutine computes the pressure and sound speed    *           
!  for all particles on a list assuming a perfect gas       *           
!  equation of state.                                       *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension rho(idim), u(idim) 
!                                                                       
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /prev/ xold(0:idim),vold(0:idim),rhold(idim),              &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /carac/ deltam(idim), abar(idim) 
      common /cgas / gamma 
!                                                                       
!--initialize quantities                                                
!                                                                       
      vsmax=0. 
      gama1=gamma - 1. 
!                                                                       
!--isothermal equation of state                                         
!                                                                       
      if(gamma.eq.1.)then 
         do k=1,ncell 
            prold(k) = pr(k) 
            pr(k)=u(k)*rho(k) 
            vsound(k)=dsqrt(pr(k)/rho(k)) 
            vsmax=dmax1(vsmax,vsound(k)) 
         enddo 
         return 
      end if 
!                                                                       
!--perfect gas equation of state                                        
!                                                                       
      do k=1,ncell 
         prold(k) = pr(k) 
         pr(k)=u(k)*gama1*rho(k) 
         vsound(k)=dsqrt(gamma*pr(k)/rho(k)) 
         vsmax=dmax1(vsound(k),vsmax) 
      enddo 
      return 
      END                                           
!                                                                       
      subroutine eospgr (ncell,rho,u) 
!************************************************************           
!                                                           *           
!  This subroutine computes the pressure, and sound speed   *           
!  according to an equation of state that includes gas and  *           
!  radiation pressure.                                      *           
!  This part of the code has not been debugged and most     *           
!  likely won't run.                                        *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension rho(idim), u(idim) 
!                                                                       
      double precision umass 
      common /prev/ xold(0:idim),vold(0:idim),rhold(idim),              &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /therm/ xmu(idim) 
      common /cgas / gamma 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
!                                                                       
      double precision rhok, uk, adrho, t9, rdmu, pgas, prad 
!                                                                       
      vsmax=0. 
      gama1=gamma - 1. 
!                                                                       
      tempmx=-1e20 
      tempmn=1e20 
      do k=1,ncell 
         rhok=dble(rho(k)) 
         uk=dble(u(k)) 
         adrho=dble(arad)/rhok 
         t9=dble(temp(k)) 
         rdmu=dble(bigr/xmu(k)) 
         call rootemp1(t9,uk,rdmu,adrho,iflag) 
         if(iflag.eq.1) then 
            write(*,*)'temperature did not converge!',k 
            write(*,*)t9,ui, rhoi 
         end if 
         temp(k)=t9 
         tempmx=dmax1(temp(k),tempmx) 
         tempmn=dmin1(temp(k),tempmn) 
!                                                                       
!--compute the various pressures                                        
!                                                                       
         prad=dble(arad)*t9*t9*t9*t9*0.3333333333333d0 
         pgas=rhoi*t9*rdmu 
         prold(k) = pr(k) 
         pr(k)=pgas+prad 
!                                                                       
!--compute thermal energy and sound speed                               
!                                                                       
         vsound(k)=dsqrt(1.66666666*real(pgas)/rho(k)) 
         vsmax=dmax1(vsound(k),vsmax) 
      enddo 
      print*,'eospgr:max,min temp(K)',tempmx*1e9,tempmn*1e9 
!                                                                       
      return 
      END        
                                         
      subroutine eos3(ncell,rho,u,ye) 
!*************************************************************          
!                                                                       
!     compute pressures and temperatures with the                       
!     Ocean eos assuming NSE                                            
!                                                                       
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      double precision umass 
      double precision rhok, uk, tempk, yek, ptot, cs, etak,            &
     &                 abark, xpk, xnk, xak, xhk, yehk, rholdk,         &
     &                 yeoldk,xpfk, p2k, p3k, p4k, xmuhk, stot          
!                                                                       
      dimension rho(idim), u(idim), ye(idim) 
!                                                                       
      common /prev/ xold(0:idim),vold(0:idim),rhold(idim),              &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /eosnu/ prnu(idim1) 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /therm/ xmu(idim) 
      common /state/ xp(idim),xn(idim),eta(idim),ifleos(idim) 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /hnucl/ xalpha(idim),xheavy(idim), yeh(idim) 
      common /cgas / gamma 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      logical ifign 
      common /ign  / ifign(idim) 
!                                                                       
      tempmx=-1e20 
      tempmn=1e20 
      vsmax=0. 
      do k=1,ncell 
         rhok=rho(k) 
         uk=u(k) 
         tempk=temp(k) 
         yek=ye(k) 
         if (yek.lt.0.) then 
            print *,'k,yek',k,yek 
            stop 
         endif 
!                                                                       
!--assume chemical freeze-out                                           
!                                                                       
         if (ifleos(k).eq.1) then 
            abark=abar(k) 
            call rootemp2(k,rhok,uk,tempk,yek,abark,                    &
     &                    ptot,cs,etak,stot)                            
            xpk=0. 
            xnk=0. 
            xmue(k)=etak*tempk 
            xmuhat(k)=0.0 
            xalpha(k)=0.0 
            xheavy(k)=1. 
            yeh(k)=yek 
            ifign(k)=.true. 
!                                                                       
!--ocean eos + nse                                                      
!                                                                       
         elseif (ifleos(k).eq.2) then 
            tempk=tempold(k) 
            xpk=xpold(k) 
            xnk=xnold(k) 
            rholdk=rhold(k)*udens 
            yeoldk=yeold(k) 
            call rootemp3(k,rhok,uk,tempk,yek,rholdk,yeoldk,            &
     &                    ptot,cs,etak,xpk,xnk,xak,xhk,yehk,abark,stot) 
            abar(k)=abark 
            xalpha(k)=xak 
            xheavy(k)=xhk 
            yeh(k)=yehk 
            xmue(k)=etak*tempk 
            xmuhat(k)=0.0 
            ifign(k)=.false.                      
!                                                                       
!--Swesty and Lattimer eos                                              
!                                                                       
         else 
            xpfk=xpf(k) 
            p2k=pvar2(k) 
            p3k=pvar3(k) 
            p4k=pvar4(k) 
            abark=abar(k) 
            call rootemp4(k,rhok,uk,yek,tempk,xpfk,p2k,p3k,p4k,         &
     &              ptot,cs,etak,xpk,xnk,xak,xhk,yehk,abark,xmuhk,stot) 
            xpf(k)=xpfk 
            pvar2(k)=p2k 
            pvar3(k)=p3k 
            pvar4(k)=p4k 
!-- assume all dissociated, this could be changed in the future         
            abar(k)=abark 
            xalpha(k)=xak 
            xheavy(k)=xhk 
            yeh(k)=yehk 
            xmue(k)=etak*tempk 
            xmuhat(k)=xmuhk 
            ifign(k)=.false. 
         endif 
!                                                                       
!--store values                                                         
!                                                                       
         if (xpk.le.1d-20) then 
            xpk=0. 
         end if 
         if (xnk.le.1d-20) then 
            xnk=0. 
         end if 
         xp(k)=xpk 
         xn(k)=xnk 
         eta(k)=etak 
         temp(k)=tempk 
         prold(k) = pr(k) 
         pr(k)=ptot 
         u2(k)=stot 
         vsound(k)=dmin1(cs,0.3333d0*clight) 
         vsmax=dmax1(vsmax,vsound(k)) 
         tempmx=dmax1(tempmx,temp(k)) 
         tempmn=dmin1(tempmn,temp(k)) 
!                                                                       
      enddo 
      !print *,'eos3: tempmn, tempmx',tempmn,tempmx                     
      !print *,'eta(1),ifleos(1),temp(1):',eta(1),ifleos(1),temp(1)     
      return 
      END                                           
!          
      
      subroutine eos5(ncell,rho,u,ye_spho) 
        !*************************************************************          
        !                                                                       
        !     compute pressures and temperatures with the                       
        !     SFHo tables                                           
        !                                                                       
        !************************************************************           
        !                                            
              use eosmodule                           
              implicit double precision (a-h,o-z) 
        !                                                                       
              parameter (idim=10000) 
              parameter (idim1=idim+1)                           
            
              real*8 xrho,xye,xtemp,xtemp2
              real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
              real*8 xdpderho,xdpdrhoe
              real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat_spho
              real*8 xxa,xxh,xxn,xxp              
              integer keytemp,keyerr
        !                                                                       
              double precision umass 
              double precision rhok, uk, tempk, yek, ptot, cs, etak,            &
             &                 abark, xpk, xnk, xak, xhk, yehk, rholdk,         &
             &                 yeoldk,xpfk, p2k, p3k, p4k, xmuhk, stot          
        !                                                                       
              dimension rho(idim), u(idim), ye_spho(idim) 
        !                                                                       
              common /prev/ xold(0:idim),vold(0:idim),rhold(idim),              &
             &              prold(idim), tempold(idim), yeold(idim),            &
             &              xnold(idim), xpold(idim)                            
              common /cases/ t9nse, rhoswe, rhonue, rhonux 
              common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
              common /eosnu/ prnu(idim1) 
              common /carac/ deltam(idim), abar(idim) 
              common /tempe/ temp(idim) 
              common /therm/ xmu(idim) 
              common /state/ xp(idim),xn(idim),eta(idim),ifleos(idim) 
              common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
              common /cpots/ xmue(idim), xmuhat(idim) 
              common /hnucl/ xalpha(idim),xheavy(idim), yeh(idim) 
              common /cgas / gamma 
              common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
              common /units/ umass, udist, udens, utime, uergg, uergcc 
              common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
              logical ifign 
              common /ign  / ifign(idim) 

              data ggcgs/6.67e-8/, avo/6.02e23/ 
              data aradcgs /7.565e-15/, boltzk/1.381e-16/ 
              data hbar/1.055e-27/, ccgs/3e10/ 
              data emssmev/0.511e0/, boltzmev/8.617e-11/ 
              data ergmev /6.2422e5/, sigma1/9d-44/, sigma2/5.6d-45/ 
              data c2cgs /6.15e-4/, c3cgs /5.04e-10/, fermi/1d-13/
        !                                                                       
              tempmx=  -1e20 !? 
              tempmn=   1e20 !?
              vsmax=    0. !?
              keytemp = 1
              keyerr  = 0
              
              urho = umass/udist**3
              upr = umass/udist/utime**2
              uv = udist/utime
              uent = 1!uergg/utemp
              ueta = ergmev/utemp
              
              print*, 'uergg ', uergg
              print*, 'urho ', urho            
              print*, 'uergg ', uergg 

              call readtable("Hempel_SFHoEOS_rho222_temp180_ye60_version_1.3_20190605.h5")

              do k=1,ncell 
                 xrho=rho(k)*urho
                 xenr=u(k)*uergg
                 xtemp=temp(k)*utemp*boltzmev
                 xye=ye_spho(k) 

                 if (xye.lt.0.) then 
                    print *,'k,yek',k,xye 
                    stop 
                 endif 
        !                                                                       
        !--SFHo EoS tables                                               
        !                     
                 call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,             &
                    xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,    &
                    xmuhat_spho,keytemp,keyerr,precision)    
                 !                                                                       
                 !--store values - CHECK UNITS!!!                                                        
                 !                 
                 abar(k)=   xabar !
                 xalpha(k)= xxa   !
                 xheavy(k)= xxh   !
                 yeh(k)=    xye   !
                 xmue(k)=   xmu_e/ergmev  !?                 
                 xmuhat(k)= xmuhat_spho/ergmev !?
                 ifign(k)=  .false. !?

                 if (xxp.le.1d-20) then 
                    xxp=0. 
                 end if 
                 if (xxn.le.1d-20) then 
                    xxn=0. 
                 end if 

                 xp(k)=     xxp !
                 xn(k)=     xxn !
                 eta(k)=    (xmu_e/xtemp)/ueta !? code units ?
                 temp(k)=   xtemp/utemp/boltzmev
                 prold(k) = pr(k) 
                 pr(k)=     xprs/upr
                 u2(k)=     xent/uent ! kb per nucleon
                 vsound(k)= sqrt(xcs2)/uv

                 ! zbar would be nice but not completely *necessary*
                 vsmax=dmax1(vsmax,vsound(k)) 
                 tempmx=dmax1(tempmx,temp(k)) 
                 tempmn=dmin1(tempmn,temp(k))                                                                     
        !                                                                       
              enddo 
              return 
              END        
                                                                        
                                                                        
      subroutine turbulence(ncell,x,f,q,v,rho,fmix) 
!****************************************************************       
!                                                               *       
!  This subroutine computes the turbulence parameters to be     *       
!  in the momentum and energy equations.    It evolves the      *       
!  turbulent velocity following Couch et al. 1902.01340         *       
!                                                               *       
!****************************************************************       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension x(0:idim),f(0:idim),v(0:idim) 
      dimension fmix(idim) 
      dimension q(idim),rho(idim) 
!                                                                       
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /fturb/ geff(idim), fmix1(idim), fmix2(idim) 
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
!                                                                       
      alpha(1)=1.d0/6.d0 
      alpha(2)=1.d0/6.d0 
      alpha(3)=1.d0/6.d0 
      alpha(4)=1.d0/6.d0 
                                                                        
      do k=1,ncell 
         dx=x(k)-x(k-1) 
         dv=v(k)-v(k-1) 
         geff(k)=geff(k)+v(k)*dv/dx 
         vturbk=dsqrt(vturb2(k)) 
         if (k.eq.1) then 
            bvf(k)=(rho(k+1)-rho(k))/dx/rho(k)-                         &
     &           (pr(k+1)-pr(k))/dx/rho(k)/vsound(k)**2                 
         else 
            bvf(k)=(rho(k+1)-rho(k))/dx/rho(k)-                         &
     &           (pr(k+1)-pr(k))/dx/rho(k)/vsound(k)**2                 
         end if 
         bvf(k)=-geff(k)*bvf(k) 
         dmix(k)=alpha(1)*pr(k)/rho(k)/geff(k) 
         dk=alpha(2)*dmix(k)*vturbk 
         if (k.eq.1) then 
            dvt=dsqrt(vturb2(k+1))-vturbk 
            dturb=1/x(k)**2/dx*(x(k)**2*rho(k+1)*                       &
     &           vturb2(k+1)*(v(k)-dmix(k+1)*alpha(2)*                  &
     &           dvt/dx) -                                              &
     &           x(k-1)**2*rho(k)*vturb2(k)*(v(k-1)-                    &
     &           dmix(k)*alpha(2)*dvt/dx))                              
         else 
            dvt=vturbk-dsqrt(vturb2(k-1)) 
            dturb=1/x(k)**2/dx*(x(k)**2*rho(k)*                         &
     &           vturb2(k)*(v(k)-dmix(k)*alpha(2)*                      &
     &           dvt/dx) -                                              &
     &           x(k-1)**2*rho(k-1)*vturb2(k-1)*(v(k-1)-                &
     &           dmix(k-1)*alpha(2)*dvt/dx))                            
         end if 
!                                                                       
!-- I'm doing a bit of a cheat here by assuming drho/dt=0.              
!                                                                       
         fmix(k)=-dturb/rho(k)-vturb2(k)*dvt+                           &
     &        vturbk*bvf(k)**2*dmix(k)-                                 &
     &        vturbk**3/dmix(k)                                         
      end do 
      return 
      END                                           
                                                                        
                                                                        
      subroutine forces(ncell,x,f,q,v,rho) 
!****************************************************************       
!                                                               *       
!  This subroutine computes the force on the cells that need to *       
!  have their forces evaluated.  We also do  neutrino diffusion *       
!  here.                                                        *       
!                                                               *       
!****************************************************************       
!                                                                       
      implicit double precision (a-h,o-z) 
!             
      logical post_bounce
      double precision pr_turb
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension x(0:idim),f(0:idim),v(0:idim) 
      dimension q(idim),rho(idim) 
!                                                                       
      common /eosnu/ prnu(idim1) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /carac/ deltam(idim), abar(idim) 
      common /damping/ damp, dcell 
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
      common /mlout/ pr_turb(idim1)
      !common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
!                                                                       
      data pi4/12.56637d0/ 
!                                                                       
!--acceleration by pressure                                             
!                                                                       
      prnu(1)=0.d0 
      prnu(ncell+1)=0.d0 
      do k=1,ncell 
         km05=k 
         kp05=k+1 
         ak=pi4*x(k)*x(k) 
         xk3=(.5d0*(x(k-1)+x(k)))**2.5 
         xkp3=(.5d0*(x(k)+x(k+1)))**2.5 
         akp1=pi4*x(k+1)*x(k+1) 
         akm1=pi4*x(k-1)*x(k-1) 
         akp05=0.5d0*(akp1+ak) 
         akm05=0.5d0*(ak+akm1) 
!         deltamk=0.5d0*(deltam(km05)+deltam(kp05))                     
!                                                                       
!--pressure gradients                                                   
!        
         pressp = pr(kp05) + prnu(kp05)          
         pressm = pr(km05) + prnu(km05)
         if (post_bounce.eqv..true.) then
             pressp = pressp+pr_turb(kp05)
             pressm = pressm+pr_turb(km05)
         endif
         gradp=ak*(pressp - pressm) 
         !gradpt=ak*(rho(kp05)*vturb2(kp05)-                             &
     !&        rho(km05)*vturb2(km05))                                   
!                                                                       
!--artificial viscosity pressure                                        
!                                                                       
         gradq=0.5d0*q(kp05)*(3.d0*akp05-ak) -                          &
     &         0.5d0*q(km05)*(3.d0*akm05-ak)                            
!                                                                       
!         write (99,199) k,x(k),f(k),                                   
!     $        gradp/deltam(km05),gradq/deltam(km05)                    
         f(k)=f(k)-(gradp+gradq+gradpt)/deltam(km05) 
!                                                                       
!--damping                                                              
!                                                                       
         if (k.lt.int(dcell)) then 
            f(k)=f(k)-damp*v(k) 
         end if 
      enddo 
! 199  format(I4,4(1x,1pe12.4))                                         
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine gravity(ncell,deltam,x,f) 
!****************************************************************       
!                                                               *       
!  This subroutine adds the gravitational force on the          *       
!  particles.  This will have the option of a neutron star      *       
!  core or it can allow the particles to make up the core.      *       
!                                                               *       
!****************************************************************       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
!                                                                       
      dimension x(0:idim),f(0:idim) 
      dimension gpot(idim),deltam(idim) 
      dimension xmi(0:idim) 
!                                                                       
      !common /fturb/ geff(idim), fmix1(idim), fmix2(idim) 
      common /core / dcore, xmcore 
      common /rshift/ gshift(idim) 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
!                                                                       
!--compute internal mass for all edges                                  
!--and calculate forces (actually accelerations)                        
!                                                                       
      if(dcore.gt.1.d0) then 
         xmi(0)=xmcore 
         r2=x(0)**2 
         f(0)=-gg*xmi(0)/r2 
      else 
         xmi(0)=0 
         f(0)=0 
      end if 
      do k=1,ncell 
         xmi(k)=xmi(k-1)+deltam(k) 
         r2=x(k)**2 
         f(k)=-gg*xmi(k)/r2 
!         geff(k)=-gg*xmi(k)/r2 
      enddo 
!                                                                       
!--calculate gravitational potential                                    
!                                                                       
      const=gg/(clight*clight) 
      r=x(ncell) 
      gpot(ncell)=-const*xmi(ncell)/r 
      rold=r 
      do k=ncell-1,1,-1 
         r=x(k) 
         gpot(k)=gpot(k+1)-(rold-r)*const*xmi(k)/(r*r) 
         rold=r 
      enddo 
!                                                                       
!--calculate gravitational redshift (w.r.t. r=infinity)                 
!                                                                       
      do k=1,ncell 
         gshift(k)=1.d0/dsqrt(1.d0-2.0d0*gpot(k)) 
!         gshift(k)=1.d0                                                
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine mmw 
!**************************************************************         
!                                                             *         
!  subroutine sets the mean molecular weight of the gas       *         
!  assuming complete ionization.                              *         
!                                                             *         
!**************************************************************         
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
!                                                                       
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /numb/ ncell, ncell1 
      common /therm/ xmu(idim) 
      common /carac/ deltam(idim), abar(idim) 
!                                                                       
      do i=1,ncell 
         xmu(i)=abar(i)/(abar(i)*ye(i)+1.) 
      enddo 
      return 
      END                                           
!                                                                       
      subroutine nserho(k,rhoold,yeold,rho,ye,t9,yp,yn,                 &
     &                  xa,xh,yeh,zbar,abar,ubind,dubind)               
!*************************************************************          
!                                                                       
! this subroutine figures out the NSE eq. assuming that yp              
! and yn were previously known at different density and ye,             
! but SAME temperature                                                  
!                                                                       
!**************************************************************         
!                                                                       
      implicit double precision(a-h,o-z) 
!                                                                       
      parameter (tolnse=1d-5,kmax=20) 
!                                                                       
!                                                                       
      common /testnse/ testk,testzy,testay,testyp,testyn 
!                                                                       
!--finding zero point of binding energy                                 
!                                                                       
      ider=2 
      call nsesolv(ider,t9,rhoold,yeold,yp,yn,kit,kmax,ubind0,          &
     &             xa,xh,yeh,zbar,abar)                                 
!                                                                       
      If(kit.ge.kmax) Then 
          write(*,*) 'NSE mis-stored entering nserho' 
          write(*,*) 'T9, rho, ye,k',t9,rhoold,yeold,k 
          write(*,*) 'inconsistent with yp, yn',yp,yn 
          stop 
      Endif 
      ypold=yp 
      ynold=yn 
      delye=ye-yeold 
      rhovar=rho-rhoold 
      ider=0 
!                                                                       
!--If delye small, skip ye variation.                                   
!                                                                       
      If(delye.lt.1.d-10) GOTO 50 
      yelast=yeold 
      Do 40 i=1,100 
          delye=dsign(min(abs(ye-yelast),abs(delye)),delye) 
          yetmp=yelast+delye 
          call nsesolv(ider,t9,rhoold,yetmp,yp,yn,kit,kmax,ubind,       &
     &                 xa,xh,yeh,zbar,abar)                             
          If (dabs(yetmp-ye).le.tolnse.and.kit.lt.kmax) goto 50 
          If (kit.ge.kmax) then 
             delye=0.5d0*delye 
             yp=ypold 
             yn=ynold 
          Elseif(kit.lt.4) Then 
             yelast=yetmp 
             delye=2.d0*delye 
             ypold=yp 
             ynold=yn 
          Else 
             yelast=yetmp 
             ypold=yp 
             ynold=yn 
          Endif 
   40 Continue 
      write(*,*)'Ye loop failure in ',i,'steps for particle',k 
      write(*,*)t9,rhoold,rho 
      write(*,*)yetmp,yelast,ye 
      write(*,*)yeold,delye,yp,yn 
      write(*,*)kit,zbar,abar 
      write(*,*)kmax,ubind,dubind 
      write(*,*)testk,testzy,testay,testyp,testyn 
      stop 
   50 Continue 
!                                                                       
!--Begin rho iteration                                                  
!                                                                       
      ypold=yp 
      ynold=yn 
      rholast=rhoold 
      If(dabs(rhovar).gt.1d7) Then 
          delrho=dsign(max(1d7,0.125d0*rhovar),rhovar) 
      Else 
          delrho=rhovar 
      Endif 
!                                                                       
      do 60 i=1,500 
          delrho=dsign(min(abs(rho-rholast),abs(delrho)),delrho) 
          rhotmp=rholast + delrho 
          call nsesolv(ider,t9,rhotmp,ye,yp,yn,kit,kmax,ubind,          &
     &                 xa,xh,yeh,zbar,abar)                             
          If(dabs((rhotmp-rho)/rho).lt.tolnse.and.kit.lt.kmax) goto 70 
          If (kit.ge.kmax) Then 
              delrho=.5d0*delrho 
              yp=ypold 
              yn=ynold 
          Elseif(kit.lt.4) Then 
             rholast=rhotmp 
             delrho=2.d0*delrho 
             ypold=yp 
             ynold=yn 
          Else 
              rholast=rhotmp 
              ypold=yp 
              ynold=yn 
          Endif 
   60 Continue 
      write(*,*)'rho loop failure in ',i,'steps for particle',k 
      write(*,*)t9,rhoold,rho,rhotmp 
      write(*,*)yetmp,yelast,ye 
      write(*,*)yeold,delye,yp,yn 
      write(*,*)kit,zbar,abar 
      write(*,*)kmax,ubind,dubind 
      write(*,*)testk,testzy,testay,testyp,testyn 
      stop 
   70 Continue 
   90 format(A25,3(1pe12.4),I3) 
!                                                                       
!--solve for actual rho and Ye, calcuate binding energy, average A and Z
!                                                                       
      ider=2 
      call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,                 &
     &             xa,xh,yeh,zbar,abar)                                 
      If(kit.ge.kmax) Then 
          write(*,*) 'Final step failure in nserho: particle ',k 
          write(*,*) kit,rhotmp,rho 
      Endif 
!                                                                       
!--Overstep in T9 to give initial estimate of dUb/dT9                   
!                                                                       
      If(t9.lt.12.d0) Then 
          delt9=1d-4 
      Elseif(t9.gt.25.d0) Then 
          delt9=1.d-1 
      Else 
          delt9=1d-2 
      Endif 
   80 t9d=t9+delt9 
      ypd=yp 
      ynd=yn 
      call nsesolv(ider,t9d,rho,ye,ypd,ynd,kit,kmax,ubindd,             &
     &             x11a,xh,yeh,zbar,abar)                               
      If(kit.ge.kmax) Then 
          write(*,*) 'dUb/dT step failure in nserho, particle ',k 
          write(*,*) kit,t9d,t9 
          delt9=.5d0*delt9 
          goto 80 
      Endif 
      dubind=(ubindd-ubind)/(t9d-t9) 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nsetemp(k,t9old,rho,ye,t9,yp,yn,                       &
     &                   xa,xh,yeh,zbar,abar,ubind,dubind)              
!*************************************************************          
!                                                                       
! this subroutine figures out the NSE eq. assuming that yp              
! and yn were previously know at the SAME density and ye,               
! but different temperature                                             
!                                                                       
!**************************************************************         
!                                                                       
!                                                                       
      implicit double precision(a-h,o-z) 
!                                                                       
      parameter (tolnse=1d-5,kmax=10) 
!                                                                       
      common /testnse/ testk,testzy,testay,testyp,testyn 
!                                                                       
!--finding zero point of binding energy                                 
!                                                                       
      ider=2 
      call nsesolv(ider,t9old,rho,ye,yp,yn,kit,kmax,ubind0,             &
     &             xa,xh,yeh,zbar,abar)                                 
!                                                                       
      If(kit.ge.kmax) Then 
          write(*,*) 'NSE mis-stored entering nsetemp' 
          write(*,*) 'T9, rho, ye',t9old,rho,ye 
          write(*,*) 'inconsistent with yp, yn',yp,yn 
          if (yp.eq.0.) stop 
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
      do i=1,1000 
          delt9=dsign(min(dabs(t9-t9last),dabs(delt9)),delt9) 
          t9tmp=t9last+delt9 
          call nsesolv(ider,t9tmp,rho,ye,yp,yn,kit,kmax,ubind,          &
     &                   xa,xh,yeh,zbar,abar)                           
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
      write(*,*)'nsetemp(1) did not converge!',k 
      write(*,*)t9old,t9tmp,t9 
      write(*,*)i,delt9,rho 
      write(*,*)ye,yp,yn 
      write(*,*)kit,kmax,ubind 
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
     &                   xa,xh,yeh,zbar,abar)                           
      If(kit.ge.kmax) Then 
          write(*,*) 'NSEtemp failed for final T9, particle',k 
          write(*,*) kit,t9,t9tmp 
          STOP 
      Endif 
!                                                                       
!--Overstep in T9 to calculate dUb/dT9                                  
!                                                                       
      If(t9.lt.12.d0) Then 
          delt9=1d-4 
      Elseif(t9.gt.25.d0) Then 
          delt9=1.d0 
      Else 
          delt9=1d-2 
      Endif 
   80 t9d=t9+delt9 
      ypd=yp 
      ynd=yn 
      call nsesolv(ider,t9d,rho,ye,ypd,ynd,kit,kmax,ubindd,             &
     &             xa,xh,yeh,zbar,abar)                                 
      If(kit.ge.kmax) Then 
          write(*,*) 'dUb/dT step failure in nsetemp, particle ',k 
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
!    &                   xa,xh,yeh,zbar,abar)                           
!     If(kit.ge.kmax) Then                                              
!         write(*,*) 'NSEtemp failed for deriv T9last: particle',k      
!         write(*,*) kit,t9,t9tmp                                       
!         STOP                                                          
!     endif                                                             
!     dubind=(ubindlast-ubind)/(t9last-t9)                              
      return 
      END                                           
!                                                                       
      subroutine nuabs(ncell,rho,x,dye,ynue,ynueb,                      &
     &                 dynue,dynueb,dunue,dunueb)                       
!****************************************************                   
!                                                                       
! this subroutine computes the neutrino absorption                      
! by nucleons.                                                          
! Note: all neutrino energies are in MeV                                
!                                                                       
!****************************************************                   
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
!                                                                       
      parameter(idim=10000) 
      parameter (delta=0.783) 
      parameter (deltab=1.805) 
!-- tffac=(6pi^2/2)^2/3 hbar*2/(2 mp kb)*avo^2/3                        
      parameter (tfermi=164.) 
!                                                                       
      dimension ctnue(0:idim),ctnueb(0:idim) 
      dimension rho(idim), dye(idim), x(0:idim) 
      dimension ynue(idim),ynueb(idim) 
      dimension dynue(idim), dynueb(idim),                              &
     &          dunue(idim),dunueb(idim)                                
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /dnuas/ dnuae(idim),dnuaeb(idim),dnuse(idim),dnuseb(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /tempe/ temp(idim) 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /rshift/ gshift(idim) 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
!                                                                       
!--the energy input is modified by delta=mn-mp-me=0.783 Mev             
!--or mp-mn-me=-1.805 Mev                                               
!                                                                       
      tol=1.d-15 
      if (enue.lt.tol) then 
         facnue=0.d0 
         xnorme=0.d0 
      else 
         facnue=(enue+delta)/enue 
         xnorme=uergg/(facnue*enue*umeverg) 
      end if 
      if (enueb.lt.tol) then 
         facnueb=0.d0 
         xnormb=0.d0 
      else 
         facnueb=(enueb-deltab)/enueb 
         xnormb=uergg/(facnueb*enueb*umeverg) 
      end if 
!                                                                       
      tf=tfermi/utemp 
      denue=0.0 
      denueb=0.0 
!      uconv=1./umevnuc                                                 
      prefac=4.*3.14159*xsecnn*clight 
      dt = steps(1) 
      dt9=9.*dt 
      do i=1,ncell 
         if (trapnue(i)) then 
!                                                                       
!--e neutrino absorption by neutrons                                    
!  ---------------------------------                                    
!                                                                       
!-- if beta equilibrium declared, set all changes to zero               
!-- particle will be taken care of in nubeta                            
            if (ebetaeq(i)) then 
               facn=0. 
            else 
!--compute degeneracy blocking                                          
               call integrals(eta(i),f0,f1,f2,f3,f4,f5,0,df2,df3) 
               call integrals(etanue(i),g0,g1,g2,g3,g4,g5,0,dg2,dg3) 
               expon=dexp(dmin1(enuet(i)/(temp(i)*utmev)-eta(i),50.d0)) 
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               tneut=dmax1(temp(i),tf*(xn(i)*rho(i)*udens)**(0.66666)) 
               tfpr=tf*(xp(i)*rho(i)*udens)**(0.66666) 
               expon=dexp((tneut-tfpr)/temp(i)) 
!--bpr: proton end-state blocking                                       
               bpr=expon/(1.+expon) 
               block=bel*bpr 
               fac=prefac*rho(i) 
               facn=fac*enuet(i)*enuet(i)*ynue(i)*xn(i)*block 
               dy9=ynue(i)/dt9 
               if (facn.gt.dy9) facn=dy9 
            endif 
            dynue(i)=dynue(i)-facn 
            enuei=enuet(i)+delta 
            facunue=facn*umevnuc*enuei 
            dnuae(i)=facunue 
            dunue(i)=dunue(i)-facunue 
            dunu(i)=dunu(i)-facunue 
            dye(i)=dye(i)+facn 
         else 
!--correct enue and e2nue for gravitational redshift                    
            fshift=1./gshift(i) 
!--artificially lower energy here                                       
            cenue=enue*fshift 
            ce2nue=e2nue*fshift*fshift 
            expon=dexp(min(cenue/(temp(i)*utmev)-eta(i),50.d0)) 
!--bel: electron end-state blocking                                     
            bel=expon/(1.+expon) 
            block=bel 
            fac=xsecnn/(x(i)*x(i)) 
!                                                                       
!  a) energy absorption                                                 
!                                                                       
            facunue=xn(i)*fac*ce2nue*rlumnue*fshift*facnue*block 
            denue=denue+facunue*deltam(i)/fshift 
            dunu(i)=dunu(i) - facunue 
!                                                                       
!  b) change in Ye                                                      
!                                                                       
            dye(i)=dye(i) + facunue*xnorme 
         endif 
!                                                                       
!--e anti-neutrino absorption by protons                                
!  -------------------------------------                                
!                                                                       
         if (trapnueb(i)) then 
            if (pbetaeq(i)) then 
               facp=0. 
            else 
!--since all the energy from the nueb goes into e+ which is never       
!--degenerate, the only term that counts is xmuhat                      
               tprot=dmax1(temp(i),tf*(xp(i)*rho(i)*udens)**(0.66666)) 
               tfne=tf*(xn(i)*rho(i)*udens)**(0.666666666) 
               expon=dexp((tprot-tfne)/temp(i)) 
!--bne: neutron end-state blocking                                      
               bne=expon/(1.+expon) 
               block=bne 
               fac=prefac*rho(i) 
               facp=fac*enuebt(i)*enuebt(i)*ynueb(i)*xp(i)*block 
               dy9=ynueb(i)/dt9 
               if (facp.gt.dy9) facp=dy9 
            endif 
            dynueb(i)=dynueb(i)-facp 
            enuebi=enuebt(i)-deltab 
            facunueb=facp*umevnuc*enuebi 
            dnuaeb(i)=facunueb 
            dunueb(i)=dunueb(i)-facunueb 
            dunu(i)=dunu(i)-facunueb 
            dye(i)=dye(i)-facp 
         else 
            fshift=1./gshift(i) 
            ce2nueb=e2nueb*fshift*fshift 
            fac=xsecnn/(x(i)*x(i)) 
!                                                                       
!  a) energy absorption                                                 
!                                                                       
            facunueb=xp(i)*fac*ce2nueb*rlumnueb*fshift*facnueb 
            denueb=denueb+facunueb*deltam(i)/fshift 
            dunu(i)=dunu(i) - facunueb 
!                                                                       
!  b) change in Ye                                                      
!                                                                       
            dye(i)=dye(i) - facunueb*xnormb 
!                                                                       
         endif 
      enddo 
!                                                                       
      fo=ufoe/utime 
      !print *, 'denue, denueb, fo', denue, denueb,fo                   
      !print*,'  e-neutrino absorption:',denue*fo,' foes/s'             
      !print*,'e-neutrino bar absorption:',denueb*fo,' foes/s'          
!                                                                       
!--check if denue larger than 0.25*rlumnue                              
!                                                                       
!--change this back to .25                                              
      if (denue.gt.0.10*rlumnue) then 
         jtrape=1 
!         print*,'nuabs: nue absorption/emission=',                     
!     1              denue/rlumnue                                      
      elseif (denue.lt.0.05*rlumnue.or.denue.lt.1.d-8) then 
         jtrape=-1 
      else 
         jtrape=0 
      endif 
!                                                                       
! --change this back to .25                                             
      if (denueb.gt.0.10*rlumnueb) then 
         jtrapb=1 
!         print*,'nuabs: nueb absorption/emission=',                    
!     1            denueb/rlumnueb                                      
      elseif (denueb.lt.0.05*rlumnueb.or.denueb.lt.1.d-8) then 
         jtrapb=-1 
      else 
         jtrapb=0 
      endif 
      return 
      END                                           
!                                                                       
      subroutine nuann(ncell,x,rho,ynue,ynueb,ynux,                     &
     &                 dynue,dynueb,dynux,dunue,dunueb,dunux)           
!************************************************************           
!                                                                       
!  This subroutine computes the rate of neutrino anti neutrino          
!  annihilation into e+/e- pairs                                        
!  see Goodman, Dar, Nussinov, ApJ 314 L7                               
!                                                                       
!*********************************************************              
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
!                                                                       
      parameter(sinw2=0.23) 
      parameter(fe=(1.+4.*sinw2+8.*sinw2*sinw2)/(6.*3.14159)) 
      parameter(fx=(1.-4.*sinw2+8.*sinw2*sinw2)/(6.*3.14159)) 
!-- cross section Gf / gram is 6.02e23*5.29e-44=3.2e-20                 
      parameter(sigmae=fe*3.2e-20) 
      parameter(sigmax=fx*3.2e-20) 
!                                                                       
      dimension rho(idim),x(0:idim) 
      dimension ynue(idim),ynueb(idim),ynux(idim) 
      dimension dynue(idim), dynueb(idim), dynux(idim),                 &
     &          dunue(idim),dunueb(idim),dunux(idim)                    
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /tempe/ temp(idim) 
      double precision umass 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
!                                                                       
      double precision dsigcue, dsigcux 
!                                                                       
      dsigcue=dble(sigmae)*umass/dble(udist*udist) 
      dsigcux=dble(sigmax)*umass/dble(udist*udist) 
      sigcue=dsigcue 
      sigcux=dsigcux 
!      uconv=1./umevnuc                                                 
      do i=1,ncell 
         tmev=utmev*temp(i) 
         dx=(x(i)-x(i-1)) 
         if (trapnux(i).and.dnux(i).lt.dx) then 
            expon=dexp(dmin1(enuxt(i)/tmev-eta(i),50.d0)) 
!--bel: electron end-state blocking                                     
            bel=expon/(1.+expon) 
            facx=sigcux*rho(i)*ynux(i)*ynux(i)*bel 
            outnux=facx*enuxt(i)*enuxt(i) 
            dynux(i)=dynux(i)-outnux 
!--4 because 4 x neutrinos                                              
            facunux=umevnuc*4.*outnux*enuxt(i) 
            dunux(i)=dunux(i)-facunux 
            dunu(i)=dunu(i)-facunux 
         endif 
!--we require BOTH nue and nueb to be trapped                           
         if (trapnue(i).and.trapnueb(i).and.                            &
     &       dnue(i).lt.dx.and.dnueb(i).lt.dx) then                     
            expon=dexp(min(0.5*(enuet(i)+enuebt(i))/tmev-               &
     &                      eta(i),50.d0))                              
!--bel: electron end-state blocking                                     
            bel=expon/(1.+expon) 
            face=sigcue*rho(i)*ynue(i)*ynueb(i)*bel 
            outnue=face*enuet(i)*enuebt(i) 
            dynue(i)=dynue(i)-outnue 
            dynueb(i)=dynueb(i)-outnue 
            facunue=umevnuc*outnue*enuet(i) 
            facunueb=umevnuc*outnue*enuebt(i) 
            dunue(i)=dunue(i)-facunue 
            dunueb(i)=dunueb(i)-facunueb 
            dunu(i)=dunu(i)-facunue-facunueb 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nubeta(ncell,x,rho,ye,dye,ynue,ynueb,unue,unueb,       &
     &                 dynue,dynueb,dunue,dunueb)                       
!*************************************************************          
!                                                                       
! This subroutine treats cases where beta eq. has occurred.             
! In beta eq.: munue(beta)=mue-muhat                                    
! so we compute Ynue(munue(beta)) and unue(munue(beta))                 
! assuming thermal distribution at matter temperature,                  
! compare with actual Ynue and unue, and move                           
! things in the right direction                                         
!                                                                       
!************************************************************           
!                                                                       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!--1/(2.*avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)               
      parameter (prefac=1.09e7) 
!--mn-mp-me in MeV                                                      
!      parameter (delta=0.783)                                          
!      parameter (delta2=delta*delta)                                   
!      parameter (deltab=1.805)                                         
!      parameter (deltab2=deltab*deltab)                                
!                                                                       
      dimension x(0:idim),rho(idim), ye(idim), dye(idim) 
      dimension ynue(idim),ynueb(idim) 
      dimension unue(idim),unueb(idim) 
      dimension dynue(idim), dynueb(idim),                              &
     &          dunue(idim),dunueb(idim)                                
!                                                                       
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /tempe/ temp(idim) 
      double precision umass 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
!                                                                       
!      umevnuct=umevnuc*utime                                           
!      uconv=1./umevnuct                                                
      yfac=prefac/udens 
      ufac=umevnuc*yfac 
      kounte=0 
      kountp=0 
      do i=1,ncell 
!-- time scale to eq. = 10 sound crossing time                          
         dx=x(i)-x(i-1) 
         dt=steps(i) 
         deltinv=0.1*dt*(vsound(i)/dx)**2 
         if (ebetaeq(i)) then 
            etabeta=eta(i)-xmuhat(i)/temp(i) 
            call integrals(etabeta,f0,f1,f2,f3,f4,f5,0,df2,df3) 
            tmev=temp(i)*utmev 
            tmev2=tmev*tmev 
            tmev3=tmev2*tmev 
            tmev4=tmev3*tmev 
            ynuebeta=yfac*f2*tmev3/rho(i) 
            unuebeta=ufac*f3*tmev4/rho(i) 
            dybeta=(ynuebeta-ynue(i))*deltinv 
            dynue(i)=dynue(i)+dybeta 
            dye(i)=dye(i)-dybeta 
            dubeta=(unuebeta-unue(i))*deltinv 
            dunue(i)=dunue(i)+dubeta 
            dunu(i)=dunu(i)+dubeta 
!            if (i.eq.1) print*,'nubeta: ye(1),dye(1),ynue(1)'          
!     $           ,ye(1),dye(1),ynue(1)                                 
         endif 
         if (pbetaeq(i)) then 
            etabeta=xmuhat(i)/temp(i)-eta(i) 
            call integrals(etabeta,f0,f1,f2,f3,f4,f5,0,df2,df3) 
            tmev=temp(i)*utmev 
            tmev2=tmev*tmev 
            tmev3=tmev2*tmev 
            tmev4=tmev3*tmev 
            ynuebeta=yfac*f2*tmev3/rho(i) 
            unuebeta=ufac*f3*tmev4/rho(i) 
            dybeta=(ynuebeta-ynueb(i))*deltinv 
            dynueb(i)=dynueb(i)+dybeta 
            dye(i)=dye(i)+dybeta 
            dubeta=(unuebeta-unueb(i))*deltinv 
            dunueb(i)=dunueb(i)+dubeta 
            dunu(i)=dunu(i)+dubeta 
!            if (i.eq.1) print*,'nubeta: ye(1),dye,ynueb(1)'            
!     $           ,ye(1),dye(1),ynueb(1)                                
            kountp=kountp+1 
         endif 
      enddo 
      !print *,'nubeta: kounte,kountp,ye(1)',kounte,kountp,ye(1)        
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nucheck(ncell,x,ye,                                    &
     &            ynue,ynueb,ynux,unue,unueb,unux)                      
!*************************************************************          
!                                                                       
!                                                                       
!  This routine determines trapping status by computing opacities       
!  and diffusion coefficients. Also mops up residual neutrinos          
!  from untrapped particles if they are present in tiny                 
!  quantities                                                           
!                                                                       
!*************************************************************          
!                                                                       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
!                                                                       
      parameter (idim=10000) 
      parameter (tiny=1d-15) 
      parameter (rcrit=1.) 
!                                                                       
      dimension x(0:idim), ye(idim) 
      dimension ynue(idim),ynueb(idim),ynux(idim),                      &
     &          unue(idim),unueb(idim),unux(idim)                       
!                                                                       
      double precision umass 
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      emop=0. 
!                                                                       
!-- check for nu trapping                                               
!                                                                       
      kountnue=0 
      kountnueb=0 
      kountnux=0 
      rnumx=0.d0 
      cmaxnue=0.0d0 
      cmaxnueb=0.0d0 
      cmaxnux=0.0d0 
      if (jtrape.eq.1) ftrape=ftrape*1.05 
      if (jtrapb.eq.1) ftrapb=ftrapb*1.05 
      if (jtrapx.eq.1) ftrapx=ftrapx*1.05 
      if (jtrape.eq.-1.and.ftrape.ge.1.05) ftrape=ftrape*0.97 
      if (jtrapb.eq.-1.and.ftrapb.ge.1.05) ftrapb=ftrapb*0.97 
      if (jtrapx.eq.-1.and.ftrapx.ge.1.05) ftrapx=ftrapx*0.97 
      !print *, 'jtrape, ftrape', jtrape, ftrape                        
!                                                                       
!      write(*,5)'nucheck: ftrape,ftrapb,ftrapx',ftrape,ftrapb,ftrapx   
    5 format(A30,3(1x,1g10.3)) 
      do k=1,ncell 
         km1=k-1 
         dx=x(k)-x(km1) 
!-- increase che-- should there be a factor here?                       
         che=ftrape*dx 
         chb=ftrapb*dx 
         chx=ftrapx*dx 
         rnue=che/dnue(k) 
         rnueb=chb/dnueb(k) 
         rnux=chx/dnux(k) 
         rnumx=dmax1(rnue,rnueb,rnux,rnumx) 
         ak=dble(k) 
!                                                                       
!--note if nux trapped, then nue trapped and contraposition             
!-rnue usually set at 1.                                                
         if (rnue.lt.rcrit) then 
            trapnue(k)=.false. 
         else 
            kountnue=kountnue+1 
            trapnue(k)=.true. 
            cmaxnue=dmax1(ak,cmaxnue) 
         endif 
         if (rnueb.lt.rcrit) then 
            trapnueb(k)=.false. 
         else 
            kountnueb=kountnueb+1 
            trapnueb(k)=.true. 
            cmaxnueb=dmax1(ak,cmaxnueb) 
         endif 
         if (rnux.lt.rcrit) then 
            trapnux(k)=.false. 
         else 
            kountnux=kountnux+1 
            trapnux(k)=.true. 
            cmaxnux=dmax1(ak,cmaxnux) 
         endif 
      enddo 
!                                                                       
!--do 2nd pass to pickup transparent particles behind the nusphere      
!                                                                       
      do k=1,ncell 
         if (.not.trapnue(k).and.k.lt.cmaxnue) then 
            trapnue(k)=.true. 
            kountnue=kountnue+1 
         endif 
         if (.not.trapnueb(k).and.k.lt.cmaxnueb) then 
            trapnueb(k)=.true. 
            kountnueb=kountnueb+1 
         endif 
         if (.not.trapnux(k).and.k.lt.cmaxnux) then 
            trapnux(k)=.true. 
            kountnux=kountnux+1 
         endif 
      enddo 
!                                                                       
!--mop up neutrino garbage                                              
!                                                                       
      do i=1,ncell 
         if (ynue(i).lt.tiny) then 
            ynue(i)=0. 
            emop=emop+unue(i)*deltam(i) 
            unue(i)=0. 
         endif 
         if (ynueb(i).lt.tiny) then 
            ynueb(i)=0. 
            emop=emop+unueb(i)*deltam(i) 
            unueb(i)=0. 
         endif 
         if (ynux(i).lt.tiny) then 
            ynux(i)=0. 
            emop=emop+unux(i)*deltam(i) 
            unux(i)=0. 
         endif 
      enddo 
!      if (emop.ne.0.0) print *,'nucheck: mopped up',emop/50.,'  foes'  
!                                                                       
!      write(*,110) rnumx, kountnue,kountnueb,kountnux                  
  110 format('nucheck: rnumx,nue trap,nueb trap,nux trap',              &
     &       1(1pe10.2), 3(1x,I4))                                      
!--total lepton number                                                  
      tlept=0.0 
      do i=1,ncell 
!         if (i.lt.10) then                                             
!            print *, i,(ye(i)+ynue(i)-ynueb(i)),trapnue(i)             
!         end if                                                        
         tlept=tlept+deltam(i)*(ye(i)+ynue(i)-ynueb(i)) 
      enddo 
      !print *,'nucheck: total lepton number:',tlept                    
!                                                                       
!--write nu emission to a file                                          
!                                                                       
!      f=ufoe/utime                                                     
!      write(59,120)t,f*rlumnue,enue,f*rlumnueb,enueb,f*rlumnux,enux    
!  120 format((1pe12.5),6(1x,1pe10.2))                                  
!                                                                       
!--write boundary position to a file                                    
!                                                                       
!      write(58,130) t                                                  
!  130 format(3(1x,1pe12.4))                                            
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nuconv(ncell,x,rho,ynue,ynueb,ynux,                    &
     &     dynue,dynueb,dynux,dunue,dunueb,dunux)                       
!**********************************************************             
!                                                                       
!     This subroutine computes the annhilation of                       
!     neutrino antineutrino pairs into neutrino antineutrino            
!     pairs of other species                                            
!     Calculation follows e+/e- scattering discussion in                
!     Mandl and Shaw, p. 315                                            
!     total cross-section is Gf^2*s/(12pi)                              
!                                                                       
!*****************************************************                  
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
      parameter(fs=1./(12.*3.14159)) 
!-- cross section Gf / gram is 6.02e23*5.29e-44=3.2e-20                 
      parameter(sigma=fs*3.2e-20) 
!                                                                       
      dimension rho(idim),x(0:idim) 
      dimension ynue(idim),ynueb(idim),ynux(idim) 
      dimension dynue(idim), dynueb(idim), dynux(idim),                 &
     &          dunue(idim),dunueb(idim),dunux(idim)                    
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /rshift/ gshift(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /tempe/ temp(idim) 
      double precision umass 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /carac/ deltam(idim), abar(idim) 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      double precision dsigcu 
!                                                                       
      dsigcu=dble(sigma)*umass/dble(udist*udist) 
      sigcu=dsigcu 
!      uconv=1./umevnuc                                                 
      do i=1,ncell 
!         tmev=temp(i)*utmev                                            
         dx=x(i)-x(i-1) 
         if (trapnux(i).and.dnux(i).lt.dx) then 
            expon=dexp(min(enuxt(i)/tempnue(i)-etanue(i),50.d0)) 
!--bnue: neutrino end-state blocking                                    
            bnue=expon/(1.+expon) 
            expon=dexp(dmin1(enuxt(i)/tempnueb(i)-etanueb(i),50.d0)) 
!--bnueb: anti-neutrino end-state blocking                              
            bnueb=expon/(1.+expon) 
            blockx=bnue*bnueb 
            sigrho=sigcu*rho(i) 
            facx=sigrho*ynux(i)*ynux(i)*blockx 
!--facnux: number of nux/nuxb annihilation into nue/nueb (per channel)  
            facnux=facx*enuxt(i)*enuxt(i) 
            outnue=-2.*facnux 
            outnux=-.5*outnue 
            dynux(i)=dynux(i)-outnux 
            facunux=facnux*umevnuc*2.*enuxt(i) 
            outunue=-facunux 
            outunueb=-facunux 
            outunux=-outunue-outunueb 
            dunux(i)=dunux(i)-outunux 
            if (facnux.ne.0.0) then 
               if (trapnue(i)) then 
!--contribution to the nue field if trapped:                            
                  dynue(i)=dynue(i)-outnue 
                  dunue(i)=dunue(i)-outunue 
               else 
!--contribution to the nue luminosities if nue not trapped:             
                  dunu(i)=dunu(i)-outunue 
                  shift=gshift(i) 
                  tle=0.5*facunux*deltam(i)*shift 
                  tenue=enuxt(i)*shift 
                  te2nue=tenue*tenue 
                  rlumnue=rlumnue+tle 
                  enue=enue+tle*tenue 
                  e2nue=e2nue+tle*te2nue 
               endif 
               if (trapnueb(i)) then 
!--contribution to the nueb field if trapped:                           
                  dynueb(i)=dynueb(i)-outnue 
                  dunueb(i)=dunueb(i)-outunueb 
               else 
!--contribution to the nueb luminosities if nueb not trapped:           
                  dunu(i)=dunu(i)-outunueb 
                  shift=gshift(i) 
                  tleb=0.5*facunux*deltam(i)*shift 
                  tenueb=enuxt(i)*shift 
                  te2nueb=tenueb*tenueb 
                  rlumnueb=rlumnueb+tleb 
                  enueb=enueb+tleb*tenueb 
                  e2nueb=e2nueb+tleb*te2nueb 
               endif 
            endif 
         endif 
!--nue-nueb conversion only if both trapped                             
         if (trapnue(i).and.trapnueb(i).and.                            &
     &       dnue(i).lt.dx.and.dnueb(i).lt.dx) then                     
            expon=dexp(min(.5*(enuet(i)+enuebt(i))/tempnux(i)-          &
     &                      etanux(i),50.d0))                           
!--bnux: xneutrino end-state blocking                                   
            bnux=expon/(1.+expon) 
            blocke=bnux*bnux 
            sigrho=sigcu*rho(i) 
            face=sigrho*ynue(i)*ynueb(i)*blocke 
!--facnue: number of nue/nueb annihilation into nux/nuxb (per channel)  
            facnue=face*enuet(i)*enuebt(i) 
            outnue=2.*facnue 
            dynue(i)=dynue(i)-outnue 
            dynueb(i)=dynueb(i)-outnue 
            facunue=facnue*umevnuc*enuet(i) 
            facunueb=facnue*umevnuc*enuebt(i) 
            outunue=facunue 
            outunueb=facunueb 
            dunue(i)=dunue(i)-outunue 
            dunueb(i)=dunueb(i)-outunueb 
!                                                                       
            if (facnue.ne.0.0) then 
               if (trapnux(i)) then 
!--contribution to the nux field if trapped:                            
                  dynux(i)=dynux(i)+0.5*outnue 
                  dunux(i)=dunux(i)+outunue 
               else 
!--contribution to the nux luminosities if nux not trapped:             
                  dunu(i)=dunu(i)-outunue-outunueb 
                  shift=gshift(i) 
                  tlx=(facunue+facunueb)*deltam(i)*shift 
                  tenux=0.5*(enuet(i)+enuebt(i))*shift 
                  te2nux=tenux*tenux 
                  rlumnux=rlumnux+tlx 
                  enux=enux+tlx*tenux 
                  e2nux=e2nux+tlx*te2nux 
               endif 
            endif 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nudiff(ncell,x,rho,time,ye,                            &
     &                 ynue,ynueb,ynux,dynue,dynueb,dynux,              &
     &                 dunue,dunueb,dunux)                              
!***************************************************************        
!                                                              *        
!  This subroutine calculates the diffusion in the trapped     *        
!  regime.  It includes a flux limiter of sorts.               *        
!                                                              *        
!***************************************************************        
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer ncell 
      parameter (tiny=1.d-10) 
      parameter (idim=10000) 
!                                                                       
      dimension x(0:idim), rho(idim), ye(idim) 
      dimension dynue(idim),dynueb(idim),dynux(idim),                   &
     &          dunue(idim),dunueb(idim),dunux(idim),                   &
     &          ynue(idim),ynueb(idim),ynux(idim)                       
!                                                                       
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /rshift/ gshift(idim) 
      common /nsat/ satc,xtime 
      common /dnuas/ dnuae(idim),dnuaeb(idim)                           &
     &     ,dnuse(idim),dnuseb(idim)                                    
      common /neutm/ iflxlm, icvb 
!                                                                       
      data pi4/12.56637d0/ 
!                                                                       
      dcnueold = 0 
      denueold = 0 
      dcnuebold = 0 
      denuebold = 0 
      dcnuxold = 0 
      denuxold = 0 
!                                                                       
!--Marc has pointed out that we shouldn't take c/3. so I have           
!--set c3 = clight.-or not                                              
!                                                                       
      c3 = clight/3. 
                                                                        
      dynuetot=0.0 
      dynuebtot=0.0 
      dynuxtot=0.0 
      do i=1,ncell 
         dynuetot=dynuetot+deltam(i)*dynue(i) 
         dynuebtot=dynuebtot+deltam(i)*dynueb(i) 
         dynuxtot=dynuxtot+deltam(i)*dynux(i) 
      enddo 
!                                                                       
!--track diffusion mechanism                                            
      knetot=0 
      knebtot=0 
      knxtot=0 
      knesat=0 
      knebsat=0 
      knxsat=0 
!                                                                       
      do i=1,ncell 
         dcnue = 0 
         denue = 0 
         dcnueb = 0 
         denueb = 0 
         dcnux = 0 
         denux = 0 
         i1 = i+1 
         dt=steps(i) 
         voli=deltam(i)/rho(i) 
         if (i1.le.ncell) then 
            voli1=deltam(i1)/rho(i1) 
            xp05=dsqrt(0.5d0*(x(i1)**2+x(i)**2)) 
            xm05=dsqrt(0.5d0*(x(i)**2+x(i-1)**2)) 
            gij=gshift(i)/gshift(i1) 
            gji=gshift(i1)/gshift(i) 
            rho1i=1.0d0/rho(i) 
            factor=umevnuc*rho1i 
            ai=pi4*x(i)*x(i) 
            rij=xp05-xm05 
            ynuei=max(ynue(i),0.d0) 
            ynuej=max(ynue(i1),0.d0) 
            ynuebi=max(ynueb(i),0.d0) 
            ynuebj=max(ynueb(i1),0.d0) 
            ynuxi=max(ynux(i),0.d0) 
            ynuxj=max(ynux(i1),0.d0) 
         else 
            voli1=voli 
            xp05=dsqrt(x(i)**2) 
            xm05=dsqrt(0.5d0*(x(i)**2+x(i-1)**2)) 
            gij=1.d0 
            gji=1.d0 
            rho1i=1.0d0/rho(i) 
            factor=umevnuc*rho1i 
            ai=pi4*x(i)*x(i) 
            rij=xp05-xm05 
            ynuei=max(ynue(i),0.d0) 
            ynuej=0.d0 
            ynuebi=max(ynueb(i),0.d0) 
            ynuebj=0.d0 
            ynuxi=max(ynux(i),0.d0) 
            ynuxj=0.d0 
            dnue(i1)=dnue(i) 
            dnueb(i1)=dnueb(i) 
            dnux(i1)=dnux(i) 
            rho(i1)=rho(i) 
            enuet(i1)=enuet(i) 
            enuebt(i1)=enuebt(i) 
            enuxt(i1)=enuxt(i) 
            tempnue(i1)=tempnue(i) 
            tempnueb(i1)=tempnueb(i) 
            tempnux(i1)=tempnux(i) 
            etanue(i1)=etanue(i) 
            etanueb(i1)=etanueb(i) 
            etanux(i1)=etanux(i) 
            trapnue(i1)=.false. 
            trapnueb(i1)=.false. 
            trapnux(i1)=.false. 
         endif 
!                                                                       
!--Now lets do each neutrino type                                       
!                                                                       
!--First, e-                                                            
!                                                                       
         if (trapnue(i)) then 
            knetot=knetot+1 
!--first take the average of the two diffusion coeffecients             
            dnueij = 2.0*(c3*dnue(i)*c3*dnue(i1))/                      &
     &          (c3*dnue(i)+c3*dnue(i1))                                
            dnueij=dnueij/rij 
            rflx=dnueij/c3 
!--numerical flux limiter                                               
            dcnueli=ynuei*rho(i)/(4.0*dt) 
            dcnuelj=ynuej*rho(i1)/(4.0*dt) 
!--effective energies                                                   
            enuei=gij*enuet(i) 
            enuej=gji*enuet(i1) 
!--blocking                                                             
            expij=dexp(dmin1(enuei/tempnue(i1) - etanue(i1),50.d0)) 
            expji=dexp(dmin1(enuej/tempnue(i) - etanue(i),50.d0)) 
            blockij = expij/(1+expij) 
            blockji = expji/(1+expji) 
            if (iflxlm.eq.1) then 
!--   wnuij is the diffusion limit (wave coefficient)                   
               wnuij = c3 
!--   Marc's "flux limiter" takes the minimum of these two              
               tnueij = dmin1(dnueij,wnuij)*ai/voli 
            else if (iflxlm.eq.2) then 
               sigr=3.d0/(1+.5d0*rflx+.125d0*rflx**2) 
               sigr=sigr+1.d0 
               dlamr=1.d0/(3.d0+rflx*sigr) 
               tnueij = 3.d0*dnueij*dlamr*ai/voli 
!               write (42,*) 'bwflx',i,tnueij,dlamr,dnueij              
            else if (iflxlm.eq.3) then 
               if (dnuae(i).lt.tiny) then 
                  dome = 1.d0 
               else 
                  dome=(dnuse(i)+blockij*dnuae(i))                      &
     &                 /(dnuse(i)+dnuae(i))                             
               end if 
               rflx=rflx/dome 
               dlamr=(2.d0+rflx)/(6.d0+3.d0*rflx+                       &
     &              rflx**2)/dome                                       
               tnueij= 3.d0*dnueij*dlamr/dome*ai/voli 
!               write(42,*) i, tnueij, dnueij*ai/voli                   
            end if 
!--effective concentrations                                             
            ceffi=rho(i)*ynuei*blockij 
            ceffj=rho(i1)*ynuej*blockji 
!--neutrino concentration changes                                       
            dcin=tnueij*ceffj 
            dcout=tnueij*ceffi 
!--numerical stability                                                  
            if (dcnuelj.lt.dcin) then 
               dcin = dcnuelj 
               knesat=knesat+1 
            end if 
            if (dcnueli.lt.dcout) then 
               dcout = dcnueli 
               knesat=knesat+1 
            end if 
!--total change in concentration and energy                             
            dcnue = dcin-dcout 
            denue = dcin*enuej - dcout*enuei 
            if (.not.trapnue(i1))then 
               if (dcnue.gt.0) then 
                  dcnue=0.d0 
                  denue=0.d0 
               else 
                  dcnue=-dcout 
                  denue=-dcout*enuei 
               endif 
            endif 
         endif 
!                                                                       
!--Now, e+                                                              
!                                                                       
         if (trapnueb(i)) then 
            knebtot=knebtot+1 
!--first take the average of the two diffusion coeffecients             
            dnuebij = 2.0*(c3*dnueb(i)*c3*dnueb(i1))/                   &
     &              (c3*dnueb(i)+c3*dnueb(i1))                          
            dnuebij = dnuebij/rij 
            rflx=dnuebij/c3 
!--numerical flux limiter                                               
            dcnuebli=ynuebi*rho(i)/(4.0*dt) 
            dcnueblj=ynuebj*rho(i1)/(4.0*dt) 
!--effective energies                                                   
            enuebi=gij*enuebt(i) 
            enuebj=gji*enuebt(i1) 
!--blocking                                                             
            expij=dexp(min(enuebi/tempnueb(i1) - etanueb(i1),50.d0)) 
            expji=dexp(min(enuebj/tempnueb(i) - etanueb(i),50.d0)) 
            blockij = expij/(1+expij) 
            blockji = expji/(1+expji) 
            if (iflxlm.eq.1) then 
!--   wnuij is the diffusion limit (wave coefficient)                   
               wnuij = c3 
!--   Marc's "flux limiter" takes the minimum of these two              
               tnuebij = dmin1(dnuebij,wnuij)*ai/voli 
!               write (42,*) 'mflx',i,tnuebij,wnuij,dnueij              
            else if (iflxlm.eq.2) then 
               sigr=3.d0/(1+.5d0*rflx+.125d0*rflx**2) 
               sigr=sigr+1.d0 
               dlamr=1.d0/(3.d0+rflx*sigr) 
               tnuebij = 3.d0*dnuebij*dlamr*ai/voli 
!               write (42,*) 'bwflx',i,tnuebij,dlamr,dnueij             
            else if (iflxlm.eq.3) then 
               if (dnuaeb(i).lt.tiny) then 
                  dome = 1.d0 
               else 
                  dome=(dnuseb(i)+blockij*dnuaeb(i))                    &
     &              /(dnuseb(i)+dnuaeb(i))                              
               end if 
               rflx=rflx/dome 
               dlamr=(2.d0+rflx)/(6.d0+3.d0*rflx+rflx**2)/dome 
               tnuebij= 3.d0*dnuebij*dlamr/dome*ai/voli 
            end if 
!--effective concentrations                                             
            ceffi=rho(i)*ynuebi*blockij 
            ceffj=rho(i1)*ynuebj*blockji 
!--neutrino concentration changes                                       
            dcin=tnuebij*ceffj 
            dcout=tnuebij*ceffi 
            dcouteb = dcout 
!--numerical stability                                                  
            if (dcnueblj.lt.dcin) then 
               dcin = dcnueblj 
               knebsat=knebsat+1 
            end if 
            if (dcnuebli.lt.dcout) then 
               dcout = dcnuebli 
               knebsat=knebsat+1 
            end if 
!--total change in concentration and energy                             
            dcnueb = dcin-dcout 
            denueb = dcin*enuebj - dcout*enuebi 
            if (.not.trapnueb(i1))then 
               if (dcnueb.gt.0) then 
                  dcnueb=0.d0 
                  denueb=0.d0 
               else 
                  dcnueb=-dcout 
                  denueb=-dcout*enuebi 
               endif 
            endif 
         endif 
!                                                                       
!--First, e-                                                            
!                                                                       
         if (trapnux(i)) then 
            knxtot=knxtot+1 
!--first take the average of the two diffusion coeffecients             
            dnuxij = 2.0*(c3*dnux(i)*c3*dnux(i1))/                      &
     &             (c3*dnux(i)+c3*dnux(i1))                             
            dnuxij = dnuxij/rij 
            rflx=dnuxij/c3 
            if (iflxlm.eq.1) then 
!--   wnuij is the diffusion limit (wave coefficient)                   
               wnuij = c3 
!--   Marc's "flux limiter" takes the minimum of these two              
               tnuxij = dmin1(dnuxij,wnuij)*ai/voli 
!               write (42,*) 'mflx',i,tnuxij,wnuij,dnueij               
            else if (iflxlm.eq.2) then 
               sigr=3.d0/(1+.5d0*rflx+.125d0*rflx**2) 
               sigr=sigr+1.d0 
               dlamr=1.d0/(3.d0+rflx*sigr) 
               tnuxij = 3.d0*dnuxij*dlamr*ai/voli 
!               write (42,*) 'bwflx',i,tnuxij,dlamr,dnueij              
            else if (iflxlm.eq.3) then 
               dome=1.d0 
               rflx=rflx/dome 
               dlamr=(2.d0+rflx)/                                       &
     &              (6.d0+3.d0*rflx+rflx**2)/dome                       
               tnuxij= 3.d0*dnuxij*dlamr/dome*ai/voli 
            end if 
!--numerical flux limiter                                               
            dcnuxli=ynuxi*rho(i)/(4.0*dt) 
            dcnuxlj=ynuxj*rho(i1)/(4.0*dt) 
!--effective energies                                                   
            enuxi=gij*enuxt(i) 
            enuxj=gji*enuxt(i1) 
!--blocking                                                             
            expij=dexp(min(enuxi/tempnux(i1) - etanux(i1),50.d0)) 
            expji=dexp(min(enuxj/tempnux(i) - etanux(i),50.d0)) 
            blockij = expij/(1+expij) 
            blockji = expji/(1+expji) 
!--effective concentrations                                             
            ceffi=rho(i)*ynuxi*blockij 
            ceffj=rho(i1)*ynuxj*blockji 
!--neutrino concentration changes                                       
            dcin=tnuxij*ceffj 
            dcout=tnuxij*ceffi 
!--numerical stability                                                  
            if (dcnuxlj.lt.dcin) then 
               dcin = dcnuxlj 
               knxsat=knxsat+1 
            end if 
            if (dcnuxli.lt.dcout) then 
               dcout = dcnuxli 
               knxsat=knxsat+1 
            end if 
!--total change in concentration and energy                             
            dcnux = dcin-dcout 
            denux = 4*(dcin*enuxj - dcout*enuxi) 
            if (.not.trapnux(i1))then 
               if (dcnux.gt.0) then 
                  dcnux=0.d0 
                  denux=0.d0 
               else 
                  dcnux=-dcout 
                  denux=-4.0d0*dcout*enuxi 
               endif 
            endif 
         endif 
!                                                                       
!--record changes... including the change due to the                    
!--inside cell.                                                         
!                                                                       
         dynuei=(dcnue-dcnueold)*rho1i 
         dynue(i)=dynue(i)+dynuei 
         dynuebi=(dcnueb-dcnuebold)*rho1i 
         dynueb(i)=dynueb(i)+dynuebi 
         dynuxi=(dcnux-dcnuxold)*rho1i 
         dynux(i)=dynux(i)+dynuxi 
!                                                                       
         dunuei=(denue-denueold)*factor 
         dunue(i)=dunue(i) + dunuei 
         dunuebi=(denueb-denuebold)*factor 
         dunueb(i)=dunueb(i) + dunuebi 
         dunuxi=(denux-denuxold)*factor 
         dunux(i)=dunux(i) + dunuxi 
!                                                                       
!-- the flux out of one cell goes into the next                         
!-- but the cells have different volumes                                
         dcnueold = dcnue*voli/voli1 
         denueold = denue*voli/voli1 
         dcnuebold = dcnueb*voli/voli1 
         denuebold = denueb*voli/voli1 
         dcnuxold = dcnux*voli/voli1 
         denuxold = denux*voli/voli1 
      enddo 
!--total lepton number                                                  
      dynuetot=0.0 
      dynuebtot=0.0 
      dynuxtot=0.0 
      do i=1,ncell 
         dynuetot=dynuetot+deltam(i)*dynue(i) 
         dynuebtot=dynuebtot+deltam(i)*dynueb(i) 
         dynuxtot=dynuxtot+deltam(i)*dynux(i) 
      enddo 
      !print *, 'total ne',knetot,'ne sat',knesat                       
      !print *, 'total neb',knebtot,'neb sat',knebsat                   
      !print *, 'total nx',knxtot,'nx sat',knxsat                       
      if (dble(knesat)/dble(max(1,knetot)).gt.0.3.or.                   &
     &     dble(knebsat)/dble(max(1,knebtot)).gt.0.3.or.                &
     &     dble(knxsat)/dble(max(1,knxtot)).gt.0.3) satc=satc+1.d0      
!      if (dble(knesat+knebsat+knxsat).lt.1) satc=satc-.5d0             
      !print *, 'isat=',satc                                            
      return 
!                                                                       
      END                                           
!                                                                       
      subroutine nuecap(ncell,rho,ye,dye,dynue,dynueb,                  &
     &                  dunue,dunueb)                                   
!****************************************************                   
!                                                                       
! this subroutine computes the neutrino production                      
! by e+/e- capture on nucleons                                          
! Note: all neutrino energies are in MeV                                
!                                                                       
!****************************************************                   
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
      parameter (idim1=idim+1) 
!--1/(avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)                  
!      parameter(prefac=2.19e7)                                         
!-- tffac=(6pi^2/2)^2/3 hbar*2/(2 mp kb)*avo^2/3                        
      parameter (tfermi=164.) 
!                                                                       
      dimension rho(idim), ye(idim), dye(idim) 
      dimension dynue(idim),dynueb(idim),dunue(idim),dunueb(idim) 
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /rshift/ gshift(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)                  
!                                                                       
      double precision umass 
!      double precision ugserg,avo                                      
!                                                                       
!      data avo/6.022d23/                                               
!                                                                       
!      ugserg=dble(utime/uergg)                                         
!      umevnuct=umevnuc*utime                                           
!      uconv=1./umevnuct                                                
!      yfac=prefac/udens                                                
      tf=tfermi/utemp 
!                                                                       
!--do all cells                                                         
!                                                                       
      dt9=9.d0*steps(1) 
      do i=1,ncell 
         tempi=temp(i) 
!         rhoi=rho(i)*udens                                             
         etai=eta(i) 
         if(tempi.ge.6..or.etai*tempi.ge.10.)then 
            deltami=deltam(i) 
            xfn=xp(i)+xn(i) 
!                                                                       
!--compute electron and positron capture                                
!  -------------------------------------                                
!                                                                       
            if(xfn.ne.0)then 
               call epcapture(erate,prate,due,dup,tempi,etai,2) 
               call integrals(eta(i),f0,f1,f2,f3,f4,f5,0,df2,df3) 
               tmev=temp(i)*utmev 
!                                                                       
!--supress rates by end state degeneracy                                
!                                                                       
               expon=dexp(dmin1(f3/f2*tmev/tempnue(i)-etanue(i),        &
     &                         50.d0))                                  
!--bnue: neutrino end-state blocking                                    
               bnue=expon/(1.+expon) 
               expon=dexp(dmin1(3.*tmev/tempnueb(i)-etanueb(i),50.d0)) 
!--bnueb: anti-neutrino end-state blocking                              
               bnueb=expon/(1.+expon) 
               tfne=tf*(xn(i)*rho(i)*udens)**(0.666666666) 
               tfpr=tf*(xp(i)*rho(i)*udens)**(0.666666666) 
               tprot=dmax1(temp(i),tfpr) 
               tneut=dmax1(temp(i),tfne) 
               expon=dexp((tprot-tfne)/temp(i)) 
!--bne: neutron end-state blocking                                      
               bne=expon/(1.+expon) 
               expon=dexp((tneut-tfpr)/temp(i)) 
!--bpr: proton end-state blocking                                       
               bpr=expon/(1.+expon) 
               blocke=bnue*bne 
               blockp=bnueb*bpr 
               erate=erate*blocke 
               prate=prate*blockp 
               due=due*blocke 
               dup=dup*blockp 
               if (prate.lt.0.) then 
                  prate=1e-10 
                  dup=0. 
               endif 
               facp=xp(i)*erate 
               facn=xn(i)*prate 
!--check for beta eq.                                                   
!               tmev2=tmev*tmev                                         
!               tmev3=tmev2*tmev                                        
!              yel=yfac*tmev3*f2/rho(i)                                 
!-- if rate*(9 time steps) is larger than the e- abundance              
!-- beta eq. is declared, nuecap zeroed-out,                            
!-- will be taken care of in nubeta                                     
!              if (erate*dt9.gt.yel.and.trapnue(i)) then                
               if (facp*dt9.gt.ye(i).and.trapnue(i)) then 
                  ebetaeq(i)=.true. 
                  facp=0.0 
                  due=0.0 
               else 
                  ebetaeq(i)=.false. 
               endif 
               call integrals(-eta(i),f0,f1,f2,f3,f4,f5,0,df2,df3) 
!               ypo=yfac*tmev3*f2/rho(i)                                
!              if (prate*dt9.gt.ypo.and.trapnueb(i)) then               
               if (facn*dt9.gt.ye(i).and.trapnueb(i)) then 
                  pbetaeq(i)=.true. 
                  facn=0.0 
                  dup=0.0 
               else 
                  pbetaeq(i)=.false. 
               endif 
!                                                                       
!--change in Ye                                                         
!                                                                       
               dye(i)=dye(i)+(facn-facp) 
!                                                                       
!--energy loss and mean neutrino energies in MeV                        
!                                                                       
               enumean=due/erate 
               enubmean=dup/prate 
               due=due*umevnuc*xp(i) 
               dup=dup*umevnuc*xn(i) 
               dunu(i)=dunu(i)+due+dup 
               if (trapnue(i)) then 
                  dynue(i)=dynue(i)+facp 
                  dunue(i)=dunue(i)+due 
!                  print *,'trap:',i,dynue(i),dunue(i)                  
               else 
!--weigh the mean energy sums by luminosity                             
                  shift=gshift(i) 
                  rlnue=due*deltami*shift 
                  rlumnue=rlumnue+rlnue 
                  enue=enue+enumean*shift*rlnue 
                  e2nue=e2nue+enumean*enumean*shift*shift*rlnue 
               end if 
               if (trapnueb(i)) then 
                  dynueb(i)=dynueb(i)+facn 
                  dunueb(i)=dunueb(i)+dup 
               else 
!--weigh the mean energy sums by luminosity                             
                  shift=gshift(i) 
                  rlnueb=dup*deltami*shift 
                  rlumnueb=rlumnueb+rlnueb 
                  enueb=enueb+enubmean*shift*rlnueb 
                  e2nueb=e2nueb+enubmean*enubmean*shift*shift*rlnueb 
               end if 
            endif 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nuinit(ncell,rho,x,ye,dye,                             &
     &            ynue,ynueb,ynux,dynue,dynueb,dynux,                   &
     &            unue,unueb,unux,dunue,dunueb,dunux)                   
!*************************************************************          
!                                                                       
!  This subroutine zeroes out whatever is necessary for the             
!  neutrino physics and computes the MeV mean energies                  
!  Note that ynux is the abundance of any single specie of              
!  tau, mu, antineutrino or neutrino, but unux is the total             
!  energy in all 4 fields                                               
!  We also compute opacities, diffusion coefficients and                
!  degeneracy for each species                                          
!                                                                       
!*************************************************************          
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
!                                                                       
!                                                                       
      parameter (idim=10000) 
      parameter (small=1d-20) 
      parameter (rcrit=1.) 
!                                                                       
!--so that we don't have to go double precision:                        
!--avo*1e-44=6.02e-20                                                   
      parameter(avosig=6.02e-21) 
!--same thing for 1e13 (fm/cm) avo**-1/3                                
      parameter(fachn=1.2e5) 
!--hbar*c/kb in fermi*Kelvin: why a cap for Kelvin and not for fermi?   
      parameter(facdeb=2.3e12) 
!--e**2/kb/(1fm)                                                        
      parameter(gamfac=1.67e10) 
!--nucleon elastic xsection from Bowers and Wilson 1982, ApJSup 50      
      parameter(sigpel=1.79) 
      parameter(signel=1.64) 
!--Mayle's thesis for alpha/nu elastic scattering                       
      parameter(sigalel=0.048) 
!--prefactor to BBAL formula (Bethe 1990) for coherent scattering       
      parameter(sigheel=1.7) 
!--neutrino absorptions by free nucleons:                               
      parameter(sigabs=9.0) 
!--neutrino-el scatterings: neutral currents                            
      parameter(sigeln=1.5) 
!--neutrino-el scatterings: charged currents                            
      parameter(sigelc=2.8) 
!--1/(avo*pi**2*(hbar*c)**3 in Mev-3 cm-3 nucleon g-1)                  
      parameter(prefac=2.19e7) 
!--factors for S3(eta)=F3(eta)+F3(-eta)                                 
      parameter(pi=3.141592) 
      parameter(s3a=7.*pi*pi*pi*pi/60.) 
      parameter(s3b=0.5*pi*pi) 
!                                                                       
      dimension rho(idim), x(0:idim) 
      dimension dynue(idim),dynueb(idim),dynux(idim),                   &
     &          dunue(idim),dunueb(idim),dunux(idim),                   &
     &          ynue(idim),ynueb(idim),ynux(idim),                      &
     &          unue(idim),unueb(idim),unux(idim)                       
      dimension ye(idim),dye(idim) 
!                                                                       
      double precision umass 
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /hnucl/ xalpha(idim),xheavy(idim), yeh(idim) 
      common /tempe/ temp(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      logical ifirst 
      save ifirst 
      data ifirst/.true./ 
!                                                                       
      data pi43/4.1887902d0/ 
!                                                                       
!--zero-out dunu, dye, dynus, denus                                     
!                                                                       
      enuemx=-1e20 
      enuebmx=-1e20 
      enuxmx=-1e20 
      uconv=1./umevnuc 
      yfac=prefac/udens 
!                                                                       
      do i=1,ncell 
         dunu(i)=0.0 
         dye(i)=0.0 
         dynue(i)=0.0 
         dynueb(i)=0.0 
         dynux(i)=0.0 
         dunue(i)=0.0 
         dunueb(i)=0.0 
         dunux(i)=0.0 
         if (ynue(i).gt.small) then 
            enuet(i)=uconv*unue(i)/ynue(i) 
            enuemx=dmax1(enuemx,enuet(i)) 
         else 
            enuet(i)=0.0 
         endif 
         if (ynueb(i).gt.small) then 
            enuebt(i)=uconv*unueb(i)/ynueb(i) 
            enuebmx=dmax1(enuebmx,enuebt(i)) 
         else 
            enuebt(i)=0.0 
         endif 
         if (ynux(i).gt.small) then 
            enuxt(i)=uconv*unux(i)/(4.*ynux(i)) 
            enuxmx=dmax1(enuxmx,enuxt(i)) 
         else 
            enuxt(i)=0.0 
         endif 
      enddo 
!                                                                       
!      write(*,115)'nuinit: max; enue,enueb,enux',                      
!     1            enuemx,enuebmx,enuxmx                                
  115 format(A30,3(1x,1pe12.4)) 
!                                                                       
      fac=avosig*udens*udist 
      do i=1,ncell 
         facrho=fac*rho(i) 
!--Debye radius                                                         
         if (eta(i)*temp(i).ge.15.) then 
            rd=10.0*facdeb/(eta(i)*utemp*temp(i)) 
         else 
            rd=89.0*dsqrt(utemp*temp(i)/(rho(i)*udens)) 
         endif 
         if (xheavy(i).ge.1e-3) then 
!--mean distance between heavy nuclei (fm)                              
            rhn=fachn/(pi43*xheavy(i)*udens*rho(i)/abar(i))**0.3333333 
!--ratio of Coulomb inter-nuclei energy to thermal energy:              
            rgamma=gamfac*(abar(i)*yeh(i))**2*dexp(-rhn/rd)/            &
     &                                  (rhn*utemp*temp(i))             
         else 
!--a big number...                                                      
            rhn=1e9 
            rgamma=0.0 
         endif 
         cohfac=abar(i)*(1.-1.08*yeh(i))/6. 
!--neutral current/heavies scattering without corrections               
         alphanh=facrho*xheavy(i)*sigheel*cohfac 
!--neutral current/nucleon and alpha  scattering                        
         alphann=facrho*(xp(i)*sigpel+xn(i)*signel+                     &
     &            xalpha(i)*sigalel)                                    
!--charged current/neutron scattering                                   
         alphacn=facrho*xn(i)*sigabs 
!--charged current/proton scattering                                    
         alphacp=facrho*xp(i)*sigabs 
!--neutral current/e+e- scattering                                      
         alphane=facrho*sigeln 
!--charged current/e+e- scattering                                      
         alphace=facrho*sigelc 
         tmev=utmev*temp(i) 
         tmev2=tmev*tmev 
!--compute local electron energies                                      
         if (eta(i)*tmev.ge.1.) then 
            etai=eta(i) 
            call integrals(etai,f0,f1,f2,f3,f4,f5,0,df2,df3) 
            elmean=tmev*f3/f2 
            elmean2=tmev2*f4/f2 
            etai2=etai*etai 
            etai4=etai2*etai2 
            yemean=yfac*tmev2*tmev2/rho(i)*(s3a+0.25*etai4+s3b*etai2) 
         else 
            elmean=3.*tmev 
            elmean2=9.*tmev2 
            yemean=ye(i)*elmean 
         endif 
         if (enuet(i).gt.1.) then 
!--figure out neutrino degeneracy                                       
            etanu=etanue(i) 
            call rooteta(i,tmev,rho(i),ynue(i),unue(i),etanu,tempnu) 
            etanue(i)=etanu 
            tempnue(i)=tempnu 
!            if (i.eq.1) print *,'ynue(1),enuet(1),etanu',              
!     1                           ynue(1),enuet(1),etanu                
            ediff=enuet(i) 
            ediff2=enuet(i)*enuet(i) 
         else 
            etanue(i)=0.0 
            tempnue(i)=tmev 
            ediff=dmax1(elmean,enue) 
            ediff2=dmax1(elmean2,e2nue) 
         endif 
!--de Broglie wavelength fm                                             
         debro=1240./ediff 
         rrd=rd/debro 
!         shield=(rrd+.5/rrd)/(rrd+1./rrd+3.)                           
         shield=1.d0 
         rrhn=rhn/debro 
!         struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)                
         struct=1.0d0 
         opacnh=struct*shield*alphanh*ediff2 
         opac=(alphann+alphacn)*ediff2+(alphane+alphace)*ediff*yemean 
         dnue(i)=1./(opac+opacnh) 
         if (enuebt(i).gt.1.) then 
            etanu=etanueb(i) 
            call rooteta(i,tmev,rho(i),ynueb(i),unueb(i),etanu,tempnu) 
            etanueb(i)=etanu 
            tempnueb(i)=tempnu 
            ediff=enuebt(i) 
            ediff2=enuebt(i)*enuebt(i) 
         else 
            etanueb(i)=0.0 
            tempnueb(i)=tmev 
            ediff=dmax1(elmean,enueb) 
            ediff2=dmax1(elmean2,e2nueb) 
         endif 
!--de Broglie wavelength fm                                             
         debro=1240./ediff 
         rrd=rd/debro 
!         shield=(rrd+0.5/rrd)/(rrd+1./rrd+3.)                          
         shield=1.0d0 
         rrhn=rhn/debro 
!         struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)                
         struct=1.0d0 
         opacnh=struct*shield*alphanh*ediff2 
         opac=(alphann+alphacp)*ediff2+(alphane+alphace)*ediff*yemean 
         dnueb(i)=1./(opac+opacnh) 
         if (enuxt(i).gt.1.) then 
            etanu=etanux(i) 
            call rooteta(i,tmev,rho(i),ynux(i),.25*unux(i),etanu,tempnu) 
            etanux(i)=etanu 
            tempnux(i)=tempnu 
            ediff=enuxt(i) 
            ediff2=enuxt(i)*enuxt(i) 
         else 
            etanux(i)=0.0 
            tempnux(i)=tmev 
            ediff=max(elmean,enux) 
            ediff2=max(elmean2,e2nux) 
         endif 
!--de Broglie wavelength fm                                             
         debro=1240./ediff 
         rrhn=rhn/debro 
!         struct=(rrhn**2+(1./rrhn)**3)/(rrhn**2+rgamma)                
         struct=1.0d0 
         opacnh=struct*alphanh*ediff2 
         opac=alphann*ediff2+alphane*ediff*yemean 
         dnux(i)=1./(opac+opacnh) 
      enddo 
!                                                                       
!--set neutrino luminosities to zero.                                   
!                                                                       
      rlumnue=0. 
      rlumnueb=0. 
      rlumnux=0. 
!                                                                       
!--initialize energy sums                                               
!                                                                       
      enue=0. 
      enueb=0. 
      enux=0. 
      e2nue=0. 
      e2nueb=0. 
      e2nux=0. 
!                                                                       
!--if first call, setup neutrino trapping flags                         
!                                                                       
      if (ifirst) then 
         ifirst=.false. 
         kountnue=0 
         kountnueb=0 
         kountnux=0 
         cmaxnue=0.0d0 
         cmaxnueb=0.0d0 
         cmaxnux=0.0d0 
         do i=1,ncell 
            dx=x(i) - x(i-1) 
            che=ftrape*dx 
            chb=ftrapb*dx 
            chx=ftrapx*dx 
            rnue=che/dnue(i) 
            rnueb=chb/dnueb(i) 
            rnux=chx/dnux(i) 
            ai=dble(i) 
!--fix rnue to 1.                                                       
            if (rnue.lt.rcrit) then 
               trapnue(i)=.false. 
            else 
               kountnue=kountnue+1 
               trapnue(i)=.true. 
               cmaxnue=dmax1(cmaxnue,ai) 
            endif 
            if (rnueb.lt.rcrit) then 
               trapnueb(i)=.false. 
            else 
               kountnueb=kountnueb+1 
               trapnueb(i)=.true. 
               cmaxnueb=dmax1(cmaxnueb,ai) 
            endif 
            if (rnux.lt.rcrit) then 
               trapnux(i)=.false. 
            else 
               kountnux=kountnux+1 
               trapnux(i)=.true. 
               cmaxnux=dmax1(cmaxnux,ai) 
            endif 
         enddo 
!                                                                       
!--do 2nd pass to pickup transparent particles behind the nusphere      
!                                                                       
         do i=1,ncell 
            if (.not.trapnue(i).and.i.lt.cmaxnue) then 
               trapnue(i)=.true. 
               kountnue=kountnue+1 
            endif 
            if (.not.trapnueb(i).and.i.lt.cmaxnueb) then 
               trapnueb(i)=.true. 
               kountnueb=kountnueb+1 
            endif 
            if (.not.trapnux(i).and.i.lt.cmaxnux) then 
               trapnux(i)=.true. 
               kountnux=kountnux+1 
            endif 
         enddo 
!         write(*,110) kountnue,kountnueb,kountnux                      
      endif 
!  110 format('nuinit: nue trapped, nueb trapped, nux trapped',         
!     1           3(1x,I4))                                             
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nulum(print_nuloss) 
!****************************************************                   
!                                                                       
! This subroutine adds the core luminosity to the                       
! the neutrinos produced by cooling processes                           
!                                                                       
!***************************************************                    
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
      logical print_nuloss 
!                                                                       
!--normalize the energy sums                                            
!                                                                       
!                                                                       
!-- artificially soften neutrinos by .5                                 
!                                                                       
      if (rlumnue.ne.0.) then 
         enue=0.8*enue/rlumnue 
         e2nue=0.64*e2nue/rlumnue 
      endif 
      if (rlumnueb.ne.0.) then 
         enueb=0.8*enueb/rlumnueb 
         e2nueb=0.64*e2nueb/rlumnueb 
      endif 
      if (rlumnux.ne.0.) then 
         enux=enux/rlumnux 
         e2nux=e2nux/rlumnux 
      endif 
!                                                                       
! divide by the mass fraction                                           
!--this should be 1.                                                    
      rlumnue=0.8*rlumnue 
      rlumnueb=0.8*rlumnueb 
      rlumnux=rlumnux 
!                                                                       
      unitf=ufoe/utime 
      if (print_nuloss.eqv..true.) then 
        write(*,510)'[nue loss (foe/s)@(MeV)]',rlumnue*unitf,           &
     &  '               ',enue
!        write(*,510)'[nueb loss(foe/s)@(MeV)]',rlumnueb*unitf,          &
!     &  '               ',enueb
!        write(*,510)'[nux loss (foe/s)@(MeV)]',rlumnux*unitf,           &
!     &  '               ',enux
  510   format(A,1p,E10.3,A,E10.3) 
      endif 
!                                                                       
      return 
      END                                             
!                                                                       
      subroutine nupp(ncell,rho,ye,dynue,dynueb,dynux,                  &
     &                dunue,dunueb,dunux)                               
!****************************************************                   
!                                                                       
! this subroutine computes the neutrino production                      
! by e+/e- capture on nucleons                                          
! Note: all neutrino energies are in MeV                                
!                                                                       
!****************************************************                   
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
!                                                                       
      dimension rho(idim), ye(idim) 
      dimension dynue(idim), dynueb(idim), dynux(idim),                 &
     &          dunue(idim),dunueb(idim),dunux(idim)                    
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /rshift/ gshift(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      data avo/6.022d23/ 
!                                                                       
      ugserg=dble(utime/uergg) 
      umevnuct=umevnuc*utime 
      uconv=1./umevnuc 
!                                                                       
!--compute pair plasma processes                                        
!  -----------------------------                                        
!                                                                       
!--loop over all particles                                              
!                                                                       
      do i=1,ncell 
         yei=ye(i) 
         tempi=temp(i) 
         rhoi=rho(i)*udens 
         if(tempi.ge.6.)then 
            deltami=deltam(i) 
            etai=eta(i) 
            rhocgs=rhoi 
            tempk=utemp*tempi 
            deta=etai 
            call pppb(rhocgs,tempk,yei,deta,dunuel,dunuxl,              &
     &                    enuel,enuxl,enuel2,enuxl2)                    
!                                                                       
!--supress rates by end state degeneracy                                
!                                                                       
            tmev=tempi*utmev 
            expon=dexp(min(enuel/tempnue(i)-etanue(i),50.d0)) 
!--bnue: e-neutrino end-state blocking                                  
            bnue=expon/(1.+expon) 
            expon=dexp(min(real(enuel)/tempnueb(i)-etanueb(i),50.d0)) 
!--bnueb: anti e-neutrino end-state blocking                            
            bnueb=expon/(1.+expon) 
            expon=dexp(min(real(enuxl)/tempnux(i)-etanux(i),50.d0)) 
!--bnux: x neutrino end-state blocking                                  
            bnux=expon/(1.+expon) 
            blocke=bnue*bnueb 
            blockx=bnux*bnux 
            dunuel=dunuel*ugserg*blocke 
            dunuxl=dunuxl*ugserg*blockx 
            dunu(i)=dunu(i)+dunuxl+dunuel 
            if (trapnue(i)) then 
               facnue=0.5*uconv*dunuel/enuel 
               dynue(i)=dynue(i)+facnue 
               dunue(i)=dunue(i)+0.5*dunuel 
            else 
               shift=gshift(i) 
               rlnue=0.5*dunuel*deltami*shift 
               rlumnue=rlumnue+rlnue 
               enue=enue+enuel*shift*rlnue 
               e2nue=e2nue+enuel2*shift*shift*rlnue 
            endif 
            if (trapnueb(i)) then 
               facnue=0.5*uconv*dunuel/enuel 
               dynueb(i)=dynueb(i)+facnue 
               dunueb(i)=dunueb(i)+0.5*dunuel 
            else 
               shift=gshift(i) 
               rlnueb=0.5*dunuel*deltami*shift 
               rlumnueb=rlumnueb+rlnueb 
               enueb=enueb+enuel*shift*rlnueb 
               e2nueb=e2nueb+enuel2*shift*shift*rlnueb 
            endif 
            if (trapnux(i)) then 
               facnux=uconv*dunuxl/enuxl 
!-- 0.25 to divide the production among 4 nu species                    
               dynux(i)=dynux(i)+0.25*facnux 
               dunux(i)=dunux(i)+dunuxl 
            else 
               shift=gshift(i) 
               rlnux=dunuxl*deltami*shift 
               rlumnux=rlumnux+rlnux 
               enux=enux+enuxl*shift*rlnux 
               e2nux=e2nux+enuxl2*shift*rlnux 
            endif 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nupress(ncell,rho,unue,unueb,unux) 
!*******************************************************                
!                                                                       
! This subroutine computes the energy trapped in the form               
! of neutrinos and derives a pressure from that                         
!                                                                       
!******************************************************                 
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension rho(idim) 
      dimension unue(idim),unueb(idim),unux(idim) 
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /eosnu/ prnu(idim1) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
!                                                                       
      ratmax=0.0 
      etanuemx=0.0 
      do i=1,ncell 
!--gamma=4/3 since nus are always relativistic                          
         rhoi3=0.3333333333333*rho(i) 
         prnui=0.0 
         if (trapnue(i)) prnui=prnui+rhoi3*unue(i) 
         if (trapnueb(i)) prnui=prnui+rhoi3*unueb(i) 
         if (trapnux(i)) prnui=prnui+rhoi3*unux(i) 
         ratmax=dmax1(prnui/pr(i),ratmax) 
         prnu(i)=prnui 
         etanuemx=dmax1(etanuemx,etanue(i)) 
      enddo 
      !print *,'nupress: maximum prnu/pr,etanuemx',ratmax,etanuemx      
!                                                                       
      return 
      END                                           
!
!                                                                       
                                                                        
      subroutine nuscat(ncell,rho,x,ynue,ynueb,ynux,                    &
     &          dunue,dunueb,dunux)                                     
!****************************************************                   
!                                                                       
! this subroutine computes the neutrino scatterings                     
! by electron and positrons.                                            
! Total cross sections from Mandl and Shaw, Quantum Field               
! Theory, p.313. I have further assumed that the mean                   
! energy transfer is 0.25*(Enu-Ee)                                      
! Note that NES is the only process which differentiates                
! populations of nux and nuxb because of the difference                 
! in electron and positron number. I have not taken that                
! into account and instead I have used fnuxm as an average              
! efficiency.                                                           
! I require the temperature to be over 6e9 K because the                
! fermi integrals are invalid below that                                
! Note: all neutrino energies are in MeV                                
!                                                                       
!****************************************************                   
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
!                                                                       
      parameter(sinw2=0.23) 
      parameter(ga=-0.5) 
      parameter(gv=2*sinw2-0.5) 
!--factor for (numu,nutau) + e- and (numub,nutaub) + e+                 
      parameter(fnux=gv*gv+ga*gv+ga*ga) 
!--factor for (numub,nutaub) + e- and (numu,nutau) + e+                 
      parameter(fnuxb=gv*gv-gv*ga+ga*ga) 
      parameter(fnuxm=0.5*(fnux+fnuxb)) 
!--factor for nue + e- and nueb + e+                                    
      parameter(fnue=(gv+1.)*(gv+1.)+(gv+1.)*(ga+1.)+(ga+1.)*(ga+1.)) 
!--factor for nueb + e- and nue + e+                                    
      parameter(fnueb=(gv+1.)*(gv+1.)-(gv+1.)*(ga+1.)+(ga+1.)*(ga+1.)) 
!                                                                       
!--a useful factor in units of /Mev**3/s                                
!--sigma(5.6e-45 cm**2)*c*(hbar**3*c**3*pi**2)**-1*(1.602e-6 ergs/Mev)**
      parameter(useful=2.2e-3) 
!                                                                       
      parameter(tiny=1d-15) 
!                                                                       
      dimension rho(idim), x(0:idim) 
      dimension ynue(idim),ynueb(idim),ynux(idim) 
      dimension dunue(idim),dunueb(idim),dunux(idim) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
      double precision umass 
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /tnus / tempnue(idim), tempnueb(idim), tempnux(idim) 
      common /rshift/ gshift(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /dnuas/ dnuae(idim),dnuaeb(idim),dnuse(idim),dnuseb(idim) 
!                                                                       
      dt=steps(1) 
      dt1=1./dt 
!      uconv=1./umevnuc                                                 
      f=ufoe/utime 
      sigfac=useful*utime 
      if (enux.ne.0) then 
         prefacx=0.5*(fnux+fnuxb)*rlumnux/enux 
      else 
         prefacx=0.0 
      endif 
      if (enue.ne.0) then 
         ratioe=rlumnue/enue 
      else 
         ratioe=0.0 
      endif 
      if (enueb.ne.0) then 
         ratioeb=rlumnueb/enueb 
      else 
         ratioeb=0.0 
      endif 
!                                                                       
      dee=0.0 
      dep=0.0 
      dex=0.0 
      do i=1,ncell 
         tempi=temp(i) 
         etai=eta(i) 
!        if (i.eq.1) print *,'nuscat: temp1,trapnue(1)',tempi,trapnue(1)
         if (tempi.ge.6.0.or.etai*tempi.ge.10.) then 
            tmev=utmev*tempi 
            tmev2=tmev*tmev 
            tmev4=tmev2*tmev2 
            temp4=tempi*tempi*tempi*tempi 
!-- electron-neutrino scattering                                        
            call integrals(etai,f0,f1,f2,f3,f4,f5,0,df2,df3) 
            elmean=tmev*f3/f2 
!--eldens is really int dne Ee                                          
            eldens=tmev4*f3 
!-- positron-neutrino scattering                                        
            call integrals(-etai,g0,g1,g2,g3,g4,g5,0,df2,df3) 
            pomean=3.*tmev 
            podens=tmev4*g2 
            fac=xsecne*temp4/(rho(i)*x(i)*x(i)) 
            deltami=deltam(i) 
!                                                                       
!--mu and tau neutrinos                                                 
!                                                                       
            if (trapnux(i).or.ynux(i).gt.tiny) then 
               expon=dexp(min((0.75*enuxt(i)+0.25*elmean)/              &
     &                          tempnux(i)-etanux(i),50.d0))            
!--benux: xneutrino end-state blocking from e- scattering               
               benux=expon/(1.+expon) 
               expon=dexp(min((0.75*elmean+0.25*enuxt(i))/              &
     &                            tmev-etai,50.d0))                     
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               expon=dexp(min((0.75*enuxt(i)+0.25*pomean)/              &
     &                            tempnux(i)-etanux(i),50.d0))          
!--bpnux: xneutrino end-state blocking from p scattering                
               bpnux=expon/(1.+expon) 
               blocke=benux*bel 
               blockp=bpnux 
               dyelnux=eldens*fnuxm*sigfac*enuxt(i)*blocke 
               qe=(1.-dexp(-dt*dyelnux))*dt1 
               dyponux=podens*fnuxm*sigfac*enuxt(i)*blockp 
               qp=(1.-dexp(-dt*dyponux))*dt1 
!-- combine both scatterings, fold into dunus                           
               facunux=umevnuc*.25*4.*ynux(i)*((enuxt(i)-elmean)*qe+    &
     &                                         (enuxt(i)-pomean)*qp)    
               dunux(i)=dunux(i)-facunux 
               dunu(i)=dunu(i)-facunux 
            endif 
            if (.not.trapnux(i)) then 
!-- heating from nux's computed from background flux                    
               shift=gshift(i) 
               fshift=1./shift 
               cenux=enux*fshift 
               ce2nux=e2nux*fshift*fshift 
!-- electron-neutrino scattering                                        
               expon=dexp(min((0.75*elmean+0.25*cenux)/tmev-            &
     &                         etai,50.d0))                             
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               blocke=bel 
               due=fac*.25*prefacx*(ce2nux*f3-cenux*tmev*f4)*blocke 
!-- positron-neutrino scattering                                        
!-- no blocking because positrons never degenerate                      
               dup=fac*.25*prefacx*(ce2nux*g3-cenux*tmev*g4) 
               dunu(i)=dunu(i)-due-dup 
               dee=dee+deltami*due*shift 
               dep=dep+deltami*dup*shift 
               dex=dex+deltami*due*shift+deltami*dup*shift 
            endif 
!                                                                       
!--electron neutrinos                                                   
!                                                                       
            if (trapnue(i)) then 
               expon=dexp(min((0.75*enuet(i)+0.25*elmean)/              &
     &                            tempnue(i)-etanue(i),50.d00))         
!--benue: eneutrino end-state blocking from e- scattering               
               benue=expon/(1.+expon) 
               expon=dexp(min((0.75*elmean+0.25*enuet(i))/tmev-         &
     &                   etai,50.d0))                                   
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               expon=dexp(min((0.75*enuet(i)+0.25*pomean)/              &
     &                          tempnue(i)-etanue(i),50.d0))            
!--bpnue: eneutrino end-state blocking from p scattering                
               bpnue=expon/(1.+expon) 
               blocke=benue*bel 
               blockp=bpnue 
!-- combine both scatterings, fold into dunus                           
               dyelnue=eldens*fnue*sigfac*enuet(i)*blocke 
               qe=(1.-dexp(-dt*dyelnue))*dt1 
               dyponue=podens*fnueb*sigfac*enuet(i)*blockp 
               qp=(1.-dexp(-dt*dyponue))*dt1 
               facunue=umevnuc*.25*ynue(i)*((enuet(i)-elmean)*qe+       &
     &                                     (enuet(i)-pomean)*qp)        
               dnuse(i)=umevnuc*.25*ynue(i)*enuet(i)*(qe+qp) 
               dunue(i)=dunue(i)-facunue 
               dunu(i)=dunu(i)-facunue 
!       if (i.eq.9) print *,'nue: benue,bel,bpnue',benue,bel,bpnue      
!       if (i.eq.9) print *,'nue: qe,qp,blocke,blockp',qe,qp,blocke,bloc
            else 
               shift=gshift(i) 
               fshift=1./shift 
               cenue=enue*fshift 
               ce2nue=e2nue*fshift*fshift 
               expon=dexp(min((0.75*elmean+0.25*cenue)/tmev-            &
     &                         etai,50.d0))                             
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               blocke=bel 
               due=fac*.25*(fnue*ratioe*(ce2nue*f3-                     &
     &                                  cenue*tmev*f4))*blocke          
               dup=fac*.25*(fnueb*ratioe*(ce2nue*g3-cenue*tmev*g4)) 
               denue=deltami*due*shift+deltami*dup*shift 
               dee=dee+deltami*due*shift 
               dep=dep+deltami*dup*shift 
               dunu(i)=dunu(i)-dup-due 
            endif 
!                                                                       
!--electron anti neutrinos                                              
!                                                                       
            if (trapnueb(i)) then 
               expon=dexp(min((0.75*enuebt(i)+0.25*elmean)/             &
     &                             tempnueb(i)-etanueb(i),50.d0))       
!--benueb: anti-eneutrino end-state blocking from e- scattering         
               benueb=expon/(1.+expon) 
               expon=dexp(min((0.75*elmean+0.25*enuebt(i))/tmev-        &
     &                   etai,50.d0))                                   
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               expon=dexp(min((0.75*enuebt(i)+0.25*pomean)/             &
     &                             tempnueb(i)-etanueb(i),50.d0))       
!--bpnueb: eneutrino end-state blocking from p scattering               
               bpnueb=expon/(1.+expon) 
               blocke=benueb*bel 
               blockp=bpnueb 
               dyelnueb=eldens*fnueb*sigfac*enuebt(i)*blocke 
               qe=(1.-dexp(-dt*dyelnueb))*dt1 
               dyponueb=podens*fnue*sigfac*enuebt(i)*blockp 
               qp=(1.-dexp(-dt*dyponueb))*dt1 
               facunueb=umevnuc*.25*ynueb(i)*((enuebt(i)-elmean)*qe+    &
     &                                        (enuebt(i)-pomean)*qp)    
               dnuseb(i)=umevnuc*.25*ynueb(i)*enuebt(i)*(qe+qp) 
               dunueb(i)=dunueb(i)-facunueb 
               dunu(i)=dunu(i)-facunueb 
!       if (i.eq.1) print *,'nueb:qe,qp,blocke,blockp',qe,qp,blocke,bloc
            else 
               shift=gshift(i) 
               fshift=1./shift 
               cenueb=enueb*fshift 
               ce2nueb=e2nueb*fshift*fshift 
               expon=dexp(min((0.75*elmean+0.25*cenueb)/tmev-           &
     &                         etai,50.d0))                             
!--bel: electron end-state blocking                                     
               bel=expon/(1.+expon) 
               blocke=bel 
               due=fac*.25*(fnueb*ratioeb*(ce2nueb*f3-                  &
     &                                     cenueb*tmev*f4))*blocke      
!-- positron neutrino scattering                                        
               dup=fac*.25*(fnue*ratioeb*(ce2nueb*g3-cenueb*tmev*g4)) 
               denueb=deltami*due*shift+deltami*dup*shift 
               dee=dee+deltami*due*shift 
               dep=dep+deltami*dup*shift 
               dunu(i)=dunu(i)-due-dup 
            endif 
         endif 
      enddo 
!                                                                       
      !print*,'e- neutrino scattering:',dee*f,' foes/s'                 
      !print*,'e+ neutrino scattering:',dep*f,' foes/s'                 
!                                                                       
!                                                                       
!--check if dex larger than 0.25*rlumnux                                
!                                                                       
      if (dex.gt.0.15*rlumnux) then 
         jtrapx=1 
!         print*,'nuscat: nux scattering=',                             
!     1              dex/rlumnux                                        
      elseif (dex.lt.0.05*rlumnux.or.dex.lt.1.d-8) then 
         jtrapx=-1 
      else 
         jtrapx=0 
      endif 
!                                                                       
      return 
      END                                           
                                                                        
      subroutine nusphere(ncell,x,                                      &
     &                  ynue,ynueb,ynux,dynue,dynueb,dynux,             &
     &                  unue,unueb,unux,dunue,dunueb,dunux)             
!************************************************************           
!                                                                       
!  This subroutine takes care of neutrino depletion from                
!  particles which are optically thin and have non-zero                 
!  ynus either because:                                                 
!     1) They have gone from optically think to optically thin          
!  or                                                                   
!     2) They have optically thick neighbours diffusing neutrinos       
!        into them                                                      
!                                                                       
!  The losses are setup so that the e-folding length is                 
!  3 time steps or the diffusion time accross the particle,             
!  whichever is larger.                                                 
!                                                                       
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=10000) 
      parameter(tiny=1d-15) 
!                                                                       
      dimension x(0:idim) 
      dimension dynue(idim),dynueb(idim),dynux(idim),                   &
     &          dunue(idim),dunueb(idim),dunux(idim),                   &
     &          ynue(idim),ynueb(idim),ynux(idim),                      &
     &          unue(idim),unueb(idim),unux(idim)                       
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /dnus/ dnue(idim),dnueb(idim),dnux(idim) 
      common /enus/ enuet(idim),enuebt(idim),enuxt(idim) 
      common /ener1/ dq(idim), dunu(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /rshift/ gshift(idim) 
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      f=ufoe/utime 
      dee=0.0 
      deeb=0.0 
      dex=0.0 
!-- minimum time scale of emission =3 time steps                        
      t3step=3.*steps(1) 
!                                                                       
      enues=0. 
      enuebs=0. 
      enuxs=0. 
      e2nues=0. 
      e2nuebs=0. 
      e2nuxs=0. 
      nkount=0 
      kount=0 
      cthird=clight/3. 
      ttranmax = 0 
      do i=1,ncell 
         dx=(x(i)-x(i-1)) 
         shift=gshift(i) 
         twave=dx/clight 
         if (.not.trapnue(i).and.ynue(i).ge.tiny) then 
            tdif=(dx)**2/(cthird*dnue(i)) 
            ttran=dmax1(tdif,twave) 
            nkount=nkount+1 
            if (ttran.lt.t3step) then 
               ttran=t3step 
               kount=kount+1 
            endif 
            ttranmax = dmax1(ttranmax,ttran) 
            fac=1./ttran 
            facye=fac*ynue(i) 
            dynue(i)=dynue(i)-facye 
            facue=fac*unue(i) 
            dunue(i)=dunue(i)-facue 
            facee=facue*deltam(i)*shift 
            rlumnue=rlumnue+facee 
            dee=dee+facee 
            enues=enues+facee*enuet(i)*shift 
            e2nues=e2nues+facee*enuet(i)*enuet(i)*shift*shift 
         endif 
         if (.not.trapnueb(i).and.ynueb(i).ge.tiny) then 
            tdif=(dx)**2/(cthird*dnueb(i)) 
            ttran=dmax1(tdif,twave) 
            if (ttran.lt.t3step) ttran=t3step 
            fac=1./ttran 
            facyeb=fac*ynueb(i) 
            dynueb(i)=dynueb(i)-facyeb 
            facueb=fac*unueb(i) 
            dunueb(i)=dunueb(i)-facueb 
            faceeb=facueb*deltam(i)*shift 
            rlumnueb=rlumnueb+faceeb 
            deeb=deeb+faceeb 
            enuebs=enuebs+faceeb*enuebt(i)*shift 
            e2nuebs=e2nuebs+faceeb*enuebt(i)*enuebt(i)*shift*shift 
         endif 
         if (.not.trapnux(i).and.ynux(i).ge.tiny) then 
            tdif=(dx)**2/(cthird*dnux(i)) 
            ttran=dmax1(tdif,twave) 
            if (ttran.lt.t3step) ttran=t3step 
            fac=1./ttran 
            facyx=fac*ynux(i) 
            dynux(i)=dynux(i)-facyx 
            facux=fac*unux(i) 
            dunux(i)=dunux(i)-facux 
            facex=facux*deltam(i)*shift 
            rlumnux=rlumnux+facex 
            dex=dex+facex 
            enuxs=enuxs+facex*enuxt(i)*shift 
            e2nuxs=e2nuxs+facex*enuxt(i)*enuxt(i)*shift*shift 
         endif 
      enddo 
      enue=enue+enues 
      e2nue=e2nue+e2nues 
      enueb=enueb+enuebs 
      e2nueb=e2nueb+e2nuebs 
      enux=enux+enuxs 
      e2nux=e2nux+e2nuxs 
      enues=enues/(dee+1e-20) 
      enuebs=enuebs/(deeb+1e-20) 
      enuxs=enuxs/(dex+1e-20) 
!                                                                       
      !print *,kount,'  fast loss cells out of',nkount                  
      !print *,'nusphere nue  losses:',dee*f,'  foes/s at',enues        
      !print *,'nusphere nueb losses:',deeb*f,'  foes/s at',enuebs      
      !print *,'nusphere nux  losses:',dex*f,'  foes/s at',enuxs        
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine nuwork(ncell,x,v,rho,                                  &
     &           unue,unueb,unux,dunue,dunueb,dunux)                    
!**********************************************************             
!                                                                       
!  This subroutine calculates the work done by the neutrino             
!  pressure and modifies the energy derivatives in consequence          
!                                                                       
!**********************************************************             
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension x(0:idim),rho(idim),v(0:idim) 
      dimension unue(idim),unueb(idim),unux(idim) 
      dimension dunue(idim),dunueb(idim),dunux(idim) 
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /eosnu / prnu(idim1) 
      common /ener1/ dq(idim), dunu(idim) 
      common /carac/ deltam(idim), abar(idim) 
!                                                                       
      data pi4/12.56637d0/ 
!                                                                       
!--compute change in specific internal energy                           
!                                                                       
      do k=1,ncell 
         ke=k-1 
         ke1=k 
         akp1=pi4*x(ke1)*x(ke1) 
         akp=pi4*x(ke)*x(ke) 
         pdv=(akp1*v(ke1)-akp*v(ke)) 
         pr3=-0.333333333*pdv*rho(k)/deltam(k) 
         if (trapnue(k)) dunue(k)=dunue(k)+unue(k)*pr3 
         if (trapnueb(k)) dunueb(k)=dunueb(k)+unueb(k)*pr3 
         if (trapnux(k)) dunux(k)=dunux(k)+unux(k)*pr3 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine pppb(rho,temp,ye,eta,dunuel,dunuxl,                    &
     &                enuel,enuxl,enuel2,enuxl2)                        
!******************************************************                 
!                                                                       
! Compute energy loss (ergs/g/s) via pair, photo, plasma and            
! bremsstahlung processes with formulae from:                           
! Itoh et al. (1989), ApJ 339, p.354                                    
! Average energies (MeV) eyeballed from:                                
! Schinder et al. (1987) , ApJ 313, p.531                               
! rho and temp should be in cgs                                         
! Note: I only implemented plasma and pair processes,                   
!       but anyone should feel free to add the rest in...               
! Note: These processes do not change Ye, neutrinos                     
!       of all families are created in pairs                            
!       denuel=denue+denueb=2*denue                                     
!       denuxl=denumu+denumub+denutau+denutaub=4*denumu                 
!                                                                       
!******************************************************                 
      implicit double precision(a-h,o-z) 
!                                                                       
      parameter(sinw2=0.23d0) 
      parameter(cv=0.5d0+2.d0*sinw2) 
      parameter(cv2=cv*cv) 
      parameter(ca=0.5d0) 
      parameter(ca2=ca*ca) 
      parameter(cvp=1.d0-cv) 
      parameter(cvp2=cvp*cvp) 
      parameter(cap=1.d0-ca) 
      parameter(cap2=cap*cap) 
      parameter(dn=2.d0) 
      parameter(facq=((cv2-ca2)+dn*(cvp2-cap2))/                        &
     &               ((cv2+ca2)+dn*(cvp2+cap2)))                        
      parameter(emassk1=1.d0/5.9302d9) 
!                                                                       
      rhoye=rho*ye 
      tme=temp*emassk1 
      tme2=tme*tme 
      tme3=tme2*tme 
      tme4=tme3*tme 
      tme6=tme4*tme2 
      tme8=tme6*tme2 
      tmep5=dsqrt(tme) 
      tmem1=1.d0/tme 
      tmem2=tmem1*tmem1 
      tmem3=tmem2*tmem1 
      xsi=(rhoye*(1.d-9))**0.3333333d0*tmem1 
      xsi2=xsi*xsi 
      xsi3=xsi2*xsi 
!--temp in Mev                                                          
      tmev=temp*8.62d-11 
!--electron chemical potential in MeV                                   
      emumev=eta*tmev 
!                                                                       
!-- pairs                                                               
!                                                                       
      if (temp.ge.1d10) then 
         fpair=dexp(-4.9924d0*xsi)*                                     &
     &         (6.002d19 + 2.084d20*xsi + 1.872d21*xsi2)/               &
     &         (xsi3 + 1.238d0*tmem1 - 0.8141d0*tmem2)                  
      else 
         fpair=dexp(-4.9924d0*xsi)*                                     &
     &         (6.002d19 + 2.084d20*xsi + 1.872d21*xsi2)/               &
     &         (xsi3 + 9.383d-1*tmem1 - 4.141d-1*tmem2 + 5.829d-2*tmem3)
      endif 
      gpair=1.d0 - 13.04d0*tme2 + 133.5d0*tme4 + 1534.d0*tme6 +         &
     &      918.6d0*tme8                                                
      qpair=((1.d0 +                                                    &
     &        rhoye/(7.692d7*tme3 + 9.715d6*tmep5))**(-0.3d0))/         &
     &      (10.748d0*tme2 + 0.3967d0*tmep5 + 1.005d0)                  
      pairfac=0.5d0*(1.d0+facq*qpair)*gpair*fpair*dexp(-2.d0*tmem1) 
      deepair=(cv2+ca2)*pairfac 
      dexpair=dn*(cvp2+cap2)*pairfac 
      epair=dmax1(6.d0*tmev,.75d0*emumev) 
!--I don't know where 6 comes from (process heavily weighted towards    
!--high energies). 0.75 because <Kinetic E of el>=3/4efermi             
      epair2=epair*epair 
!                                                                       
!-- plasma                                                              
!                                                                       
      fplas=dexp(-0.56457d0*xsi)*                                       &
     &      (2.32d-7 + 8.449d-8*xsi + 1.787d-8*xsi2)/                   &
     &      (xsi3 + 2.581d-2*tmem1 + 1.734d-2*tmem2 + 6.99d-4*tmem3)    
      plasfac=rhoye*rhoye*rhoye*fplas 
      deeplas=cv2*plasfac 
      dexplas=dn*cvp2*plasfac 
!--why 0.05d0, don't know have to trust Schindler                       
!--tmev because that's the energy of the photons                        
      eplas=dmax1(0.5d0*tmev,0.05d0*emumev) 
      eplas2=eplas*eplas 
!                                                                       
!--total emission ergs/cc/s                                             
!                                                                       
      denuel=deepair+deeplas 
      denuxl=dexpair+dexplas 
!                                                                       
!--emission ergs/g/s                                                    
!                                                                       
      dunuel=denuel/rho 
      dunuxl=denuxl/rho 
!                                                                       
!--luminosity weighted energies:                                        
!                                                                       
      enuel=(deepair*epair+deeplas*eplas)/denuel 
      enuel2=(deepair*epair2+deeplas*eplas2)/denuel 
      enuxl=(dexpair*epair+dexplas*eplas)/denuxl 
      enuxl2=(dexpair*epair2+dexplas*eplas2)/denuxl 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine preset 
!****************************************************************       
!                                                               *       
!  This subroutine sets up all quantities needed before         *       
!  starting a simulation.                                       *       
!                                                               *       
!****************************************************************       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (idim=10000) 
!                                                                       
!                                                                       
      common /carac/ deltam(idim), abar(idim) 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /numb/ ncell,ncell1 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common /therm/ xmu(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /uswest/ usltemp, uslrho, uslu, uslp, u2sluncell1 
      common /nsat/ satc,xtime 
      common /heat/ iheat 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!--set up mean molecular weight according to matter types               
!                                                                       
      call mmw 
!                                                                       
!--set the critical densities for Swesty's eos, nu trapping             
!                                                                       
      t9nse=8. 
!      rhoswe=1.e9/udens 
      rhoswe=1.e11/udens                                               
      rhonue=0.1e11/udens 
      rhonux=0.2e11/udens 
!                                                                       
      satc=0.d0 
      xtime=5.d-5 
!                                                                       
!--remove heated cells                                                  
      iheat=7 
!                                                                       
!--get nuclear data                                                     
!--nucdata is a subroutine in nse5.f                                    
!                                                                       
      call nucdata 
!                                                                       
!--initialize the swesty eos (loadmx is in sleos.f)                     
!--note that this assumes that the density are known for the dump       
!                                                                       
      call loadmx 
      do 5 i=1,ncell 
         xpf(i)=rho(i)*uslrho 
         pvar2(i)=0.155 
         pvar3(i)=-15. 
         pvar4(i)=-10. 
         etanue(i)=0. 
         etanueb(i)=0. 
         etanux(i)=0. 
    5 continue 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine rmp(t,ncell,x,v,u,rho,ye,rb,                           &
     &         ynue,ynueb,ynux,unue,unueb,unux,nrem)                    
!************************************************************           
!                                                           *           
!  subroutine to remove particles that cross the neutron    *           
!  star's surface. Their mass is added to the neutron star  *           
!  mass.                                                    *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
      parameter(idim=10000) 
!                                                                       
      dimension x(0:idim), v(0:idim),u(idim),                           &
     &     rho(idim), ye(idim)                                          
      dimension ynue(idim),ynueb(idim),ynux(idim),                      &
     &          unue(idim),unueb(idim),unux(idim)                       
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /prev / xold(0:idim),vold(0:idim),rhold(idim),             &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /core / dcore, xmcore 
      common /heat/ iheat 
!                                                                       
      nrem=0 
!                                                                       
      rhocrit=3.d13 
      tol=1.d-15 
      tmacc=0. 
      if (rho(1).gt.rhocrit.and.ye(1).lt.0.1)then 
         tmacc=tmacc+deltam(1) 
         ncell=ncell-1 
         nrem=1 
         rb=x(1) 
         do i=1,ncell 
            x(i)=x(i+1) 
            v(i)=v(i+1) 
            deltam(i)=deltam(i+1) 
            u(i)=u(i+1) 
            xold(i)=xold(i+1) 
            rhold(i)=rhold(i+1) 
            tempold(i)=tempold(i+1) 
            yeold(i)=yeold(i+1) 
            xnold(i)=xnold(i+1) 
            xpold(i)=xpold(i+1) 
            rho(i)=rho(i+1) 
            ye(i)=ye(i+1) 
            xn(i)=xn(i+1) 
            xp(i)=xp(i+1) 
            ifleos(i)=ifleos(i+1) 
            trapnue(i)=trapnue(i+1) 
            trapnueb(i)=trapnueb(i+1) 
            trapnux(i)=trapnux(i+1) 
            eta(i)=eta(i+1) 
            temp(i)=temp(i+1) 
            abar(i)=abar(i+1) 
            istep(i)=istep(i+1) 
            steps(i)=steps(i+1) 
            t0(i)=t0(i+1) 
            if (i.eq.1) then 
               if (ynue(i).gt.tol) then 
                  unue(i)=(ynue(i)*unue(i)+                             &
     &                 ynue(i+1)*unue(i+1))/(ynue(i)+ynue(i+1))         
                  ynue(i)=ynue(i)+ynue(i+1) 
               else 
                  ynue(i)=ynue(i+1) 
                  unue(i)=unue(i+1) 
               end if 
               if (ynueb(i).gt.tol) then 
                  unueb(i)=(ynueb(i)*unueb(i)+                          &
     &                 ynueb(i+1)*unueb(i+1))/(ynueb(i)+ynueb(i+1))     
                  ynueb(i)=ynueb(i)+ynueb(i+1) 
               else 
                  ynueb(i)=ynueb(i+1) 
                  unueb(i)=unueb(i+1) 
               end if 
               if (ynux(i).gt.tol) then 
                  unux(i)=(ynux(i)*unux(i)+                             &
     &                 ynux(i+1)*unux(i+1))/(ynux(i)+ynux(i+1))         
                  ynux(i)=ynux(i)+ynux(i+1) 
               else 
                  ynux(i)=ynux(i+1) 
                  unux(i)=unux(i+1) 
               end if 
            else 
               ynue(i)=ynue(i+1) 
               ynueb(i)=ynueb(i+1) 
               ynux(i)=ynux(i+1) 
               unue(i)=unue(i+1) 
               unueb(i)=unueb(i+1) 
               unux(i)=unux(i+1) 
            end if 
         enddo 
      end if 
!                                                                       
!--add accreted mass to the neutron star                                
!                                                                       
      xmcore=xmcore+tmacc 
      return 
      END                                           
!                                                                       
      subroutine amp(t,ncell,x,v,u,rho,ye,rb,                           &
     &         ynue,ynueb,ynux,unue,unueb,unux)                         
!************************************************************           
!                                                           *           
!  subroutine to divide cells if their sizes become too     *           
!  large.  Currently, nothing is done to equilibrate        *           
!  forces.                                                  *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
      parameter(idim=10000) 
      parameter(idim1=idim+1) 
!                                                                       
      dimension x(0:idim), v(0:idim),u(idim),                           &
     &     rho(idim), ye(idim)                                          
      dimension ynue(idim),ynueb(idim),ynux(idim),                      &
     &          unue(idim),unueb(idim),unux(idim)                       
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /prev / xold(0:idim),vold(0:idim),rhold(idim),             &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /carac/ deltam(idim), abar(idim) 
      common /tempe/ temp(idim) 
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /core / dcore, xmcore 
      common /heat/ iheat 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
!                                                                       
      pi43=4.18879 
!                                                                       
      nrem=0 
!                                                                       
      xcrit=9.d-4 
      xvcrit=5.d-1 
      ncr=8 
      xvalf=0.d0 
!      do j=ncr,ncr+3                                                   
      do j=ncr,ncell 
         xval=1.d0-x(j-1)/x(j) 
         if (xval.gt.xvalf) then 
            jf=j 
            xvalf=xval 
         end if 
      end do 
      if (xvalf.gt.xvcrit) then 
         ncell=ncell+1 
         pr(ncell+1)=pr(ncell) 
         idiv=jf 
         xj=.5d0*(x(jf-1)+x(jf)) 
         do i=ncell,idiv+1,-1 
            v(i)=v(i-1) 
            xold(i)=xold(i-1) 
            rhold(i)=rhold(i-1) 
            yeold(i)=yeold(i-1) 
            xnold(i)=xnold(i-1) 
            xpold(i)=xpold(i-1) 
            rho(i)=rho(i-1) 
            x(i)=x(i-1) 
            if (i.eq.idiv+1) then 
               deltam(i)=pi43*rho(i)*(x(i)**3-xj**3) 
               tempold(i)=.9*tempold(i-1) 
               temp(i)=.95*temp(i-1) 
               u(i)=.95*u(i-1) 
            else 
               deltam(i)=deltam(i-1) 
               tempold(i)=tempold(i-1) 
               temp(i)=temp(i-1) 
               u(i)=u(i-1) 
            end if 
            ye(i)=ye(i-1) 
            xn(i)=xn(i-1) 
            xp(i)=xp(i-1) 
            ifleos(i)=ifleos(i-1) 
            trapnue(i)=trapnue(i-1) 
            trapnueb(i)=trapnueb(i-1) 
            trapnux(i)=trapnux(i-1) 
            eta(i)=eta(i-1) 
            abar(i)=abar(i-1) 
            istep(i)=istep(i-1) 
            steps(i)=steps(i-1) 
            t0(i)=t0(i-1) 
            ynue(i)=ynue(i-1) 
            ynueb(i)=ynueb(i-1) 
            ynux(i)=ynux(i-1) 
            unue(i)=unue(i-1) 
            unueb(i)=unueb(i-1) 
            unux(i)=unux(i-1) 
         enddo 
         x(jf)=xj 
         deltam(jf)=pi43*rho(jf)*(x(jf)**3-x(jf-1)**3) 
      end if 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine rootemp1(temp,utr,rdmu,adrho,iflag) 
!***************************************************************        
!                                                              *        
!     subroutine computes temperature using a Newton-Raphson   *        
!     procedure found in the numerical recipes, p.254          *        
!                                                              *        
!***************************************************************        
      implicit double precision (a-h,o-z) 
!                                                                       
      data itmax/80/ , tol/1.d-2/ 
!                                                                       
!--find root allowing a maximum of itmax iteration                      
!                                                                       
      iflag=0 
      do i=1,itmax 
         at93=adrho*temp*temp*temp 
         df=4.d0*at93 + 1.5d0*rdmu 
         f=at93*temp + 1.5d0*rdmu*temp - utr 
         dtemp=f/df 
         temp=temp - dtemp 
         if(abs(dtemp/temp).lt.tol) return 
      enddo 
      iflag=1 
!                                                                       
!--iteration went out of bound or did not converge                      
!                                                                       
      print *,'rootemp: iteration did not converge'                    
      write(*,*)temp,dtemp 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine rootemp2(i,rhoi,ui,tempi,yei,                          &
     &                    abar,ptot,cs,eta,stot)                        
!*****************************************************************      
!                                                                *      
!  Given rho, u, ye and an initial T, this                       *      
!  subroutine iterates over temperatures with a Newton-Raphson   *      
!  scheme coupled with bissection to determine thermodynamical   *      
!  variables.                                                    *      
!  The ocean eos is called to determine perfect gas, radiation   *      
!  and electron/positron contributions, assuming complete        *      
!  ionization. The coulomb subroutine determines Coulomb         *      
!  corrections                                                   *      
!                                                                *      
!*****************************************************************      
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (itmax=100) 
      parameter (dtol=1.d-2) 
!                                                                       
!      real*4 t9nse, rhoswe, rhonue, rhonux                             
!                                                                       
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
!                                                                       
!--change to the appropriate units                                      
!                                                                       
!      write (40,*) i,rhoi,ui,tempi,yei,abar                            
      t9=tempi*uotemp1 
      rho=rhoi*uorho1 
      u=ui*uou1 
!      write (41,*) t9,uotemp1,rho,uorho1,u,uou1                        
!                                                                       
!--set brackets for T iterations                                        
!                                                                       
      t9l=1d-4 
      t9h=3d1 
!                                                                       
!--compute zbar from abar and ye                                        
!                                                                       
      zbar=yei*abar 
!                                                                       
!--compute coulomb correction (since coulomb corr. not dependant        
!  on T, and freeze-out is assumed call only once)                      
!                                                                       
      call coulomb(rhoi,zbar,yei,ucoul,pcoul) 
      ucoul=ucoul*uou1 
!                                                                       
!--use Newton-Raphson to find T                                         
!                                                                       
      call nados(t9,rho,zbar,abar,pel,eel,sel,                          &
     &           ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)               
      ures=ucoul+etot-u 
      dut=det 
!                                                                       
      dt9=0.d0 
      do k=1,itmax 
         dt9old=dt9 
         if (((t9-t9h)*dut-ures)*                                       &
     &      ((t9-t9l)*dut-ures).ge.0.d0.or.                             &
     &      dabs(2.d0*ures).gt.dabs(dt9old*dut)) then                   
            dt9=0.5d0*(t9h-t9l) 
            t9=t9l+dt9 
         else 
            dt9=ures/dut 
            t9=t9-dt9 
         endif 
         if (dabs(dt9/t9).lt.dtol) goto 20 
         call nados(t9,rho,zbar,abar,pel,eel,sel,                       &
     &              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)            
         ures=etot+ucoul-u 
         dut=det 
         if (ures.lt.0.d0) then 
            t9l=t9 
         else 
            t9h=t9 
         endif 
      enddo 
!                                                                       
!--did not converge, print out error message and stop                   
!                                                                       
!      print *,'rootemp2: no convergence for part. i,rho(cgs)',         
!     1         i,rhoi                                                  
      stop 
!                                                                       
!--iteration sucessful, transform back in code units                    
!                                                                       
   20 continue 
      tempi=t9*uotemp 
      ptot=ptot*uopr+pcoul 
      cs=dsqrt(gamm*ptot/rhoi) 
      stot=stot*uotemp1/uou1 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine rootemp3(k,rhok,uk,tempk,yek,rhoold,yeold,             &
     &                    ptot,cs,eta,yp,yn,xa,xh,yeh,abar,stot)        
!*****************************************************************      
!                                                                *      
!  Given rho, u, an old T, rho, yp and yn this                   *      
!  subroutine iterates over temperatures with a Newton-Raphson   *      
!  scheme coupled with bissection to determine thermodynamical   *      
!  variables.                                                    *      
!  The ocean eos is called to determine perfect gas, radiation   *      
!  and electron/positron contributions, assuming complete        *      
!  ionization. The coulomb subroutine determines Coulomb         *      
!  corrections. The nserho and nsetemp subroutines compute       *      
!  abar, zbar, and the amount of energy gone into dissociating   *      
!  nuclei.                                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (itmax=80) 
      parameter (dtol=1d-2) 
!                                                                       
!      real  udist, udens, utime, uergg, uergcc                         
!                                                                       
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
!                                                                       
!--change to the appropriate units                                      
!                                                                       
      t9=tempk*uotemp1 
      t9old=t9 
      rho=rhok*uorho1 
      u=uk*uou1 
      t9l=1d-1 
      t9h=2d2 
!                                                                       
!--transform in cgs for the NSE table                                   
!                                                                       
      rhocgs=rhok*udens 
!      ucgs=uk*uergg                                                    
      call nserho(k,rhoold,yeold,rhocgs,yek,t9,yp,yn,                   &
     &            xa,xh,yeh,zbar,abar,ubind,dubind)                     
!      xfn=yp+yn                                                        
!                                                                       
!--transform back in Ocean's units                                      
!                                                                       
      ediss=ubind/uergg*uou1 
      dediss=dubind/uergg*uou1 
!                                                                       
!--add coulomb correction                                               
!                                                                       
      call coulomb(rhok,zbar,yek,ucoul,pcoul) 
      ucoul=ucoul*uou1 
!                                                                       
!--use Newton-Raphson to iterate for T                                  
!                                                                       
!--note abar is averaged over massive nuclei, zbar overall              
      abar2=zbar/yek 
      call nados(t9,rho,zbar,abar2,pel,eel,sel,                         &
     &           ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)               
      ures=(etot+ediss+ucoul)-u 
      dut=det+dediss 
      dt9old=dabs(t9h-t9l) 
      dt9=dt9old 
      do iter=1,itmax 
         if (((t9-t9h)*dut-ures)*                                       &
     &      ((t9-t9l)*dut-ures).ge.0.d0.or.                             &
     &      dabs(2.d0*ures).gt.dabs(dt9old*dut)) then                   
            dt9old=dt9 
            t9old=t9 
            dt9=0.5d0*(t9h-t9l) 
            t9=t9l+dt9 
         else 
            dt9old=dt9 
            t9old=t9 
            dt9=ures/dut 
            t9=t9-dt9 
!-- this is because coming down from high temperatures,                 
!-- dediss can be too small                                             
            if (t9.lt.0.5d0*t9old) t9=0.5d0*t9old 
         endif 
         if (t9.lt.2.d0) then 
            print *,'k,t9,t9old',k,t9,t9old 
            print *,'u,ures,dut',u,ures,dut 
            print *,'rhocgs,rhoold',rhocgs,rhoold 
            print *,'yek,yeold',yek,yeold 
            print*,'ediss,etot,ucoul',ediss,etot,ucoul 
            stop 
         endif 
!                                                                       
!--solve nse                                                            
!                                                                       
         call nsetemp(k,t9old,rhocgs,yek,t9,yp,yn,                      &
     &                xa,xh,yeh,zbar,abar,ubind,dubind)                 
!                                                                       
!--if tolerance is satisfied, exit                                      
!                                                                       
         if (dabs(dt9/t9).lt.dtol) goto 20 
!                                                                       
!--get derivatives for new iteration                                    
!                                                                       
         call coulomb(rhok,zbar,yek,ucoul,pcoul) 
         ucoul=ucoul*uou1 
         abar2=zbar/yek 
         call nados(t9,rho,zbar,abar2,pel,eel,sel,                      &
     &              ptot,etot,stot,dpt,det,dpd,ded,gamm,eta)            
!                                                                       
         ediss=ubind/uergg*uou1 
         dediss=dubind/uergg*uou1 
         ures=(etot+ediss+ucoul)-u 
         dut=det+dediss 
         if (ures.lt.0.d0) then 
            t9l=t9 
         else 
            t9h=t9 
         endif 
      enddo 
!                                                                       
!--did not converge, print out error message and stop                   
!                                                                       
      print *,'rootemp3: no convergence for part. k,rho(cgs)',          &
     &         k,rhok                                                   
      stop 
!                                                                       
!--iteration sucessful, transform back in code units                    
!                                                                       
   20 continue 
      tempk=t9*uotemp 
      ptot=ptot*uopr+pcoul 
      cs=dsqrt(gamm*ptot/rhok) 
      stot=stot*uotemp1/uou1 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine rootemp4(k,rhok,uk,yek,tempk,xpfk,p2k,p3k,p4k,         &
     &                   press,cs,eta,yp,yn,xa,xh,yeh,abar,xmuhk,stot)  
!**************************************************************         
!                                                                       
!     This subroutine computes the temperature                          
!     with Doug Swesty's eos using                                      
!     a Newton-Raphson procedure coupled with                           
!     bissection to prevent convergence problems.                       
!                                                                       
!**************************************************************                 
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(itmax=80) 
      parameter(dtol=1d-2) 
!                                                                       
!      real utemp, utmev, ufoe, umevnuc, umeverg                        
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu 
!                                                                       
      double precision inpvar(4) 
!                                                                       
      templ=.1d0 
      temph=1d3 
      dtemp=dabs(temph-templ) 
!-- all the u's in this routine will be MeV/baryons                     
!-- all the u's in this routine will be kB/baryons                      
      u=uk*uslu 
      temp=tempk*usltemp 
      inpvar(1)=temp 
      inpvar(2)=p2k 
      inpvar(3)=p3k 
      inpvar(4)=p4k 
      pprev=xpfk 
      brydns=rhok*uslrho 
      call slwrap(k,inpvar,yek,brydns,pprev,                            &
     &      psl,usl,dusl,gamsl,eta,yp,yn,xa,xh,yeh,abar,xmuh,u2sl)      
      ures=usl-u 
      do 10 i=1,itmax 
         dtempold=dtemp 
         if (((temp-temph)*dusl-ures)*                                  &
     &      ((temp-templ)*dusl-ures).ge.0.d0.or.                        &
     &      dabs(2.d0*ures).gt.dabs(dtempold*dusl)) then                
            dtemp=0.5d0*(temph-templ) 
            temp=templ+dtemp 
         else 
            dtemp=ures/dusl 
            temp=temp-dtemp 
         endif 
         if (dabs(dtemp/temp).lt.dtol) goto 20 
         inpvar(1)=temp 
         call slwrap(k,inpvar,yek,brydns,pprev,                         &
     &            psl,usl,dusl,gamsl,eta,yp,yn,xa,xh,yeh,abar,xmuh,u2sl)
         ures=usl-u 
         if (ures.lt.0.0d0) then 
            templ=temp 
         else 
            temph=temp 
         endif 
   10 continue 
      print *,'rootemp4: no convergence for particle k, rho',           &
     &         k,rhok                                                   
   20 continue 
      p2k=inpvar(2) 
      p3k=inpvar(3) 
      p4k=inpvar(4) 
      xpfk=pprev 
!-- convert back to code units                                          
      press=psl/uslp 
      tempk=temp/usltemp 
      cs=dsqrt(max(gamsl,1.d0)*press/rhok) 
!--xmuhat is eta*T with T in code units                                 
      xmuhk=xmuh/utmev 
      stot=u2sl/u2slu 
!      if(k.eq.1)write(*,100)'root4(1):gamma,s,xh,yeh,abar',            
!     1                            gamsl,usl,xh,yeh,abar                
  100 format(A28,5(1g10.3)) 
      return 
      END                                           
!                    
!      
      subroutine rooteta(k,tmev,rho,ynu,unu,eta,temp) 
!***********************************************************            
!                                                                       
!  subroutine computes the degree of degeneracy in a nu-field           
!  assuming a thermalized Fermi distribution, i.e.:                     
!     n = ynu*rho ~ F2(eta)                                             
!     E = unu*rho ~ F3(eta)                                             
!                                                                       
!***********************************************************            
!                                                                       
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
! const= 2*pi**2*hbar**3*c**3                                           
!--2.*avo*pi**2*(hbar*c)**3 in Mev+3 cm+3 nucleon g-1)                  
      parameter(prefac=9.2e-8) 
      parameter(etamin=1e-2) 
      parameter(etamax=50.0) 
      parameter(tol=1e-3) 
!                                                                       
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
!                                                                       
      yfac=prefac*udens 
      ufac=yfac/umevnuc 
      a=(rho*ynu*yfac)**(0.333333333333) 
      b=dsqrt(dsqrt(rho*unu*ufac)) 
      if (eta.lt.etamin) eta=etamin 
!                                                                       
      do i=1,80 
         call integrals(eta,f0,f1,f2,f3,f4,f5,1,df2,df3) 
         f3_14=dsqrt(dsqrt(f3)) 
         f2_13=f2**(0.33333333333) 
!         t3=a*f3_14                                                    
!         t2=b*f2_13                                                    
         f=a*f3_14 - b*f2_13 
         df=a*0.25*f3_14/f3*df3 - b*0.33333333*f2_13/f2*df2 
         deta=f/df 
         eta=eta-deta 
!--"effective" neutrino temperature in MeVs                             
         temp=0.5*(a/f2_13+b/f3_14) 
         if (eta.ge.etamax) then 
            eta=etamax 
            return 
         elseif (eta.lt.etamin) then 
            eta=0. 
            temp=tmev 
            return 
         elseif (abs(deta/eta).lt.tol) then 
            return 
         endif 
      enddo 
!                                                                       
      print *,'eta convergence failed: rho, ynu, unu, k',               &
     &        rho, ynu, unu, k                                          
      stop 
!                                                                       
      END                                           
!                                                                       
      subroutine slwrap(kcell,inpvar,yesl,brydns,pprev,                 &
     &     psl,usl,dusl,gamsl,etasl,                                    &
     &     ypsl,ynsl,xasl,xhsl,yehsl,abar,xmuh,u2sl)                    
!******************************************************************     
!                                                                       
!  This is the wrapper routine for the swesty-lattimer eos.             
!                                                                       
!******************************************************************     
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
!                                                                       
!-- 0.46*(mn-mp-me) + bindfe56 MeV/nucleon                              
!     parameter(ushift=0.46d0*0.783d0+8.7904d0)                         
      parameter(ushift=0.0d0) 
!                                                                       
      common /typef/ iextf, ieos 
!                                                                       
      include 'eos_m4a.inc' 
      include 'el_eos.inc' 
!                                                                       
!--these variables are needed but not declared in the include files     
      integer sf 
      double precision told, pprev, u2sl 
      double precision yesl, psl, usl, dusl, gamsl, etasl, ypsl, ynsl 
      double precision xasl, xhsl, yehsl, abar, xmuh 
!                                                                       
      ye=dmax1(yesl,0.031d0) 
      if (ye.gt..5d0) print *, 'kcell, ye: ',kcell,ye 
      ye=dmin1(ye,0.5d0) 
      call inveos(inpvar,told,ye,brydns,1,eosflg,0,sf,                  &
     &            xprev,pprev)                                          
!                                                                       
      if (sf.ne.1) print *,'inveos fails for particle',kcell 
      psl=ptot 
      if (ieos.eq.3) then 
         usl=utot+ushift 
         dusl=dudt 
         u2sl=stot 
      else 
!--entropy as variable of state                                         
         usl=stot 
         dusl=dsdt 
         u2sl=utot+ushift 
      endif 
!                                                                       
      gamsl=gam_s 
      etasl=musube/inpvar(1) 
      xmuh=muhat 
!                                                                       
!-- free (exterior) nucleon fractions                                   
      ypsl=xprot 
      ynsl=xnut 
!--mass fraction and proton fraction of heavy nuclei                    
      xhsl=xh 
      yehsl=x 
      xasl=xalfa 
      abar=a 
!                                                                       
      return 
      END
!
!                                                                       
      subroutine readini
!*************************************************************          
!                                                                       
!     reads initial conditions                                          
!                                                                       
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx,mlin_grid_size,idump,skip_dump
      logical from_dump, post_bounce      
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (iqn=17) 
      real ycc,yccave
      character*128 mlmodel_name
!                                                                       
      common /cc   / ycc(idim,iqn), yccave(iqn) 
      common /ener1/ dq(idim), dunu(idim) 
      common /numb/ ncell, ncell1 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux     
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /eosq/ pr(idim1), vsound(idim), u2(idim), vsmax 
      common /freez/ ufreez(idim) 
      common /carac/ deltam(idim), abar(idim) 
      common /core / dcore, xmcore 
      common /typef/ iextf, ieos 
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /shock/ cq,cl 
      common /tempe/ temp(idim) 
      common /cgas/ gamma 
      common /timej/ time, dt
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x
      common /propt/ dtime,tmax 
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /ener2/ tkin, tterm 
      common /outp/ rout, p1out, p2out 
      logical te(idim), teb(idim), tx(idim) 
      equivalence(trapnue,te) 
      equivalence(trapnueb,teb) 
      equivalence(trapnux,tx) 
      common /nuprint/ nups,nupk,tacr 
      common /damping/ damp,dcell 
      common /neutm/ iflxlm, icvb 
      common /ufactor/ ufact,yefact 
      logical ifign
      common /ign  / ifign(idim) 
      !common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim)
      common /mlmod/ mlmodel_name
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
      &               dum2v(idim)
      common /cent/ dj(idim)
      common /dump/ from_dump
      common /idump/ idump
      common /interp/ mlin_grid_size
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
!                                                                       
      character*1024 filin,filout 

      data pi4/12.56637d0/ 
      gg=13.34 
      tacr=1.d2 
      from_dump=.false.
!                                                                       
!--read options                                                         
!                                                                       
      open(11,file='setup') 
      read(11,*)
      read(11,10) filin 
      read(11,*)
      read(11,10) filout 
      read(11,*)
      read(11,10) mlmodel_name
   10 format(A) 
      read(11,*)
      read(11,*) mlin_grid_size
      read(11,*)      
      read(11,*) idump 
      read(11,*)
      read(11,*) dtime
      read(11,*)
      read(11,*) tmax
      read(11,*)
      read(11,*) cq
      read(11,*)
      read(11,*) cl 
      read(11,*)
      read(11,*) iextf
      read(11,*)
      read(11,*) ieos
      read(11,*)
      read(11,*) dcore 
      read(11,*)
      read(11,*) delp
      read(11,*)
      read(11,*) nups
      read(11,*)
      read(11,*) damp
      read(11,*)
      read(11,*) dcell 
      read(11,*)
      read(11,*) iflxlm
      read(11,*)
      read(11,*) icvb
      read(11,*)
      read(11,*) ufact
      read(11,*)
      read(11,*) yefact 
      !print *,'cq,cl',cq,cl 
      !print *,'iextf,ieos',iextf,ieos 
      print*,'================ Setup ================'
      print*, 'Input File:            ', trim(filin)
      print*, 'Output File:           ', trim(filout)
      print*, 'ML Model Name:         ', trim(mlmodel_name)
      print*, 'Grid Size for ML:      ', mlin_grid_size
      print*, 'Dump # to read:        ', idump
      print*, 'Dump time interval (s):', dtime
      print*, 'Max time (s):          ', tmax
      !print*, cq,cl
      !print*, iextf,ieos,dcore
      !print*, ncell,delp,nups,damp,dcell
      !print*, iflxlm, icvb, ufact, yefact
!                                                                       
!--convert to code time units                      
!        
      dtime = dtime/10
      tmax = tmax/10 
!                                                                       
!--open binary file containing initial conditions                       
!                                                             
      open(60,file=trim(filin),form='unformatted') 
      open(61,file=trim(filout),form='unformatted')
!
!--adjust position pointer relative to individual binary file
!             
      read(60) idummy        
      skip_dump=idump-idummy
      close(60)
!                                                                       
!--position pointer in binary file                                      
!               
      open(60,file=trim(filin),form='unformatted')        
      do i=1,skip_dump
         read(60) idummy                  
      enddo 
!                                                                       
!--read data                                                            
!                                                                       
      nqn=17 
!       
      read(60) idump,nc,t,xmcore,rb,ftrape,ftrapb,ftrapx,              &
            pns_ind,pns_x,shock_ind,shock_x,                           &
            bounce_time,from_dump,rlumnue,rlumnueb,rlumnux,           &
            (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),  &
            (u(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),         &
            (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),           &
            (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),          &
            (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),       &
            (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),       &
            (ufreez(i),i=1,nc),(pr(i),i=1,nc),(u2(i),i=1,nc),          &
            (dj(i),i=1,nc),                                            &
            (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),             &
            (steps(i),i=1,nc),((ycc(i,j),j=1,nqn),i=1,nc)
!        (vturb2(i),i=1,nc),                                          &
!
      print*, 'idump as read = ', idump
      
      if (shock_ind.ne.0) then
            post_bounce = .true.
      else
            post_bounce = .false.
      endif

      ncell = nc
      time = t 
!                                                                       
      !print *, nc, t, xmcore, rb, ftrape,ftrapb,ftrapx 
      !do k=1,nc 
      !   write(44,*)k,temp(k),u(k),u2(k) 
      !end do 
      do k=1,nc 
         if (ifleos(k).lt.0.6) stop 
         if (ifleos(k).eq.3) then 
            ifign(k)=.false. 
         elseif (ifleos(k).eq.2) then 
            ifign(k)=.false. 
         end if 
!         u(k)=0.9*u(k)                                                 
!         u2(k)=0.9*u2(k)                                               
      end do 
      !do k=1,3 
      !   print *,u2(k),pr(k),ufreez(k) 
      !   print *, unue(k),unueb(k),unux(k) 
      !   print *, ynue(k),ynueb(k),ynux(k) 
      !   print *, xp(k), xn(k), ifleos(k), xmcore 
      !   do j=1,nqn 
      !      print *,ycc(k,j) 
      !   end do 
      !end do 
      p1out=pr(ncell) 
      p2out=1.39*gg*deltam(ncell)/x(ncell)**4/pi4 
      rout=x(ncell) 
      
      print*,'================ ',trim(filin),' ================='
      print *, 'ncell = ', ncell 
      print *, 'xmcore =     ', xmcore 
      print *, 'rho(1) =     ', rho(1) 
      print *, 'x(ncell) =   ',x(ncell)
      print *, 'pr(ncell) =  ',pr(ncell) 
      print *, 'deltam(1) =  ', deltam(1), ' (mass of the cell, cell centered)' 
      print *, 'time =       ',time
      print*,'======================================='
!      gamma=1.666666666666667                                          
      ncell1=ncell+1 
      if (ieos.eq.4) then 
         do i=1,ncell 
            !write(71,*)i,u2(i),u(i),abar(i) 
            if (ifleos(i).eq.3) then 
               s=u2(i) 
               u2(i)=u(i) 
               u(i)=s 
               !write (71,*) i,u2(i),u(i),abar(i) 
            endif 
         enddo 
      endif 
!                                                                       
      do i=1,ncell 
         if (unueb(i).lt.0.) then 
             unueb(i)=0. 
             ynueb(i)=0.0 
         endif 
         if (unux(i).lt.0.) then 
             unux(i)=0. 
             ynux(i)=0.0 
         endif 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine printout(lu) 
!**************************************************************         
!                                                             *         
!  This subroutine prints out all the results.                *         
!                                                             *         
!**************************************************************         
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
      logical from_dump
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (iqn=17) 
!                                                                       
      real ycc,yccave 
!                                                                       
      common /cc   / ycc(idim,iqn), yccave(iqn) 
      common /ener1/ dq(idim), dunu(idim) 
      common /numb/ ncell, ncell1
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /freez/ ufreez(idim) 
      common /core / dcore, xmcore 
      common /typef/ iextf, ieos 
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /tempe/ temp(idim) 
      common /cgas / gamma 
      common /timej / time, dt 
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /carac/ deltam(idim), abar(idim) 
      common /ener2/ tkin, tterm
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)
      common /cent/ dj(idim)
      common /dump/ from_dump
      common /idump/ idump
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
      !common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
      logical te(idim), teb(idim), tx(idim) 
      dimension uint(idim), s(idim) 
      equivalence(trapnue,te) 
      equivalence(trapnueb,teb) 
      equivalence(trapnux,tx) 
      
      from_dump=.true.
!                                                                       
      if (ieos.eq.4) then 
         do i=1,ncell 
            if (ifleos(i).eq.3) then 
               s(i)=u(i) 
               uint(i)=u2(i) 
            else 
               s(i)=u2(i) 
               uint(i)=u(i) 
            endif 
         enddo 
      else 
         do i=1,ncell 
            uint(i)=u(i) 
         enddo 
      endif 
      t=time 
      nc=ncell 
      nqn=17                                                                                                                            
!       
      write(lu) idump,nc,t,xmcore,rb,ftrape,ftrapb,ftrapx,             &
            pns_ind,pns_x,shock_ind,shock_x,                           &
            bounce_time,from_dump,rlumnue,rlumnueb,rlumnux,            &
            (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),  &
            (uint(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),      &
            (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),           &
            (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),          &
            (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),       &
            (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),       &
            (ufreez(i),i=1,nc),(pr(i),i=1,nc),(s(i),i=1,nc),           &
            (dj(i),i=1,nc),                                            &
            (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),             &
            (steps(i),i=1,nc),((ycc(i,j),j=1,nqn),i=1,nc)             
!          (vturb2(i),i=1,nc),                                          &
      !print *, nc,t,xmcore,rb,ftrape,ftrapb,ftrapx                     
!               
      return 
!                                                                       
!--an error as occured while writting                                   
!                                                                       
   10  print *,'an error has occured while writing' 
!                                                                       
                                                                        
      return 
      END                                           
      subroutine integrals(eta,f0,f1,f2,f3,f4,f5,ider,df2,df3) 
!************************************************************           
!                                                           *           
!  This subroutine calculates numerical approximations to   *           
!  the fermi integrals (Takahashi, El Eid, Hillebrandt, 1978*           
!  AA, 67, 185.                                             *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
!                                                                       
!                                                                       
!--case where eta > 1e-3                                                
!                                                                       
      if(abs(eta).gt.1.d-3) then 
         eta2=eta*eta 
         eta3=eta*eta2 
         eta4=eta*eta3 
         eta5=eta*eta4 
         eta6=eta*eta5 
         if (eta.le.30.) then 
            f1=(0.5*eta2+1.6449)/                                       &
     &         (1.+dexp(-1.6855*eta))                                   
            f2num=eta3/3.+3.2899*eta 
            f2den=1.-dexp(-1.8246*eta) 
            f2=f2num/f2den 
            f3num=0.25*eta4+4.9348*eta2+11.3644 
            f3den=1.+dexp(-1.9039*eta) 
            f3=f3num/f3den 
            f4=(0.2*eta5+6.5797*eta3+45.4576*eta)/                      &
     &         (1.-dexp(-1.9484*eta))                                   
            f5=(eta6/6.+8.2247*eta4+113.6439*eta2+236.5323)/            &
     &         (1.+dexp(-1.9727*eta))                                   
            f0=eta+dlog(1.+dexp(-eta)) 
            if (ider.eq.1) then 
               df2=((eta2+3.2899)*f2den -                               &
     &              f2num*1.8246*dexp(-1.8246*eta))/                    &
     &             (f2den*f2den)                                        
               df3=((eta3+9.8696*eta)*f3den +                           &
     &              f3num*1.9039*dexp(-1.9039*eta))/                    &
     &             (f3den*f3den)                                        
            endif 
         else 
            f1=0.5*eta2+1.6449 
            f2=eta3/3.+3.2899*eta 
            f3=0.25*eta4+4.9348*eta2+11.3644 
            f4=0.2*eta5+6.5797*eta3+45.4576*eta 
            f5=eta6/6.+8.2247*eta4+113.6439*eta2+236.5323 
            f0=eta 
            if (ider.eq.1) then 
               df2=eta2+3.2899 
               df3=eta3+9.8696*eta 
            endif 
         endif 
      elseif (eta.lt.-40.) then 
         f0=eta 
         f1=0.0 
         f2=0.0 
         f3=0.0 
         f4=0.0 
         f5=0.0 
      else 
         expeta=dexp(eta) 
         f1=expeta/(1.+0.2159*dexp(0.8857*eta)) 
         f2=2.*expeta/(1.+0.1092*dexp(0.8908*eta)) 
         f3=6.*expeta/(1.+0.0559*dexp(0.9069*eta)) 
         f4=24.*expeta/(1.+0.0287*dexp(0.9257*eta)) 
         f5=120.*expeta/(1.+0.0147*dexp(0.9431*eta)) 
         f0=eta+dlog(1.+dexp(-eta)) 
      end if 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine epcapture(erate,prate,due,dup,tempi,eta,iflag) 
!************************************************************           
!                                                           *           
!  subroutine to compute the energy loss due to electron    *           
!  capture on protons. Formalism taken from Takahashi, El   *           
!  Eid, Hillebrandt, 1978, AA, 67, 185. (The energy release *           
!  is in Mev/utime/nucleon and the rate in /utime/nucleon   *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter (delta=1.531) 
      parameter (q1=2.-2.*delta) 
      parameter (q2=(1.-8.*delta+2.*delta**2)/2.) 
      parameter (q3=2.*delta**2-delta) 
      parameter (q4=(4.*delta**2-1.)/8.) 
      parameter (p1=(2.-3.*delta)) 
      parameter (p2=(1.-12.*delta+6.*delta**2)/2.) 
      parameter (p3=(-3.*delta+12.*delta**2-2.*delta**3)/2.) 
      parameter (p4=(-1.+12.*delta**2-16.*delta**3)/8.) 
      parameter (p5=(2.+3*delta**2-4.*delta**3)/8.) 
      common /epcap/ betafac, c2cu, c3cu 
      common /neutm/ iflxlm, icvb 
!                                                                       
      beta=betafac/tempi 
      beta2=beta*beta 
      beta3=beta*beta2 
      beta4=beta*beta3 
      beta5=beta*beta4 
      beta6=beta*beta5 
!                                                                       
!--compute fermi integrals for electrons                                
!                                                                       
      call integrals(eta,fe0,fe1,fe2,fe3,fe4,fe5,0,df2,df3) 
!                                                                       
!--compute fermi integrals for positrons                                
!                                                                       
      call integrals(-eta,fp0,fp1,fp2,fp3,fp4,fp5,0,df2,df3) 
!                                                                       
!--e- + p -> n + nu rate / nucleons                                     
!                                                                       
      if (icvb.eq.1) then 
         erate=fe4 + q1*fe3*beta + q2*fe2*beta2 + q3*fe1*beta3 +        &
     &        q4*fe0*beta4                                              
      else 
!-- Cooperstein, Van den Horn, Baron method                             
         erate=fe4 
      end if 
      erate=c2cu*erate/beta5 
!                                                                       
!--e+ + n -> p + nubar rate / nucleon                                   
!                                                                       
      if (icvb.eq.1) then 
         prate=fp4 - q1*fp3*beta + q2*fp2*beta2 - q3*fp1*beta3 +        &
     &        q4*fp0*beta4                                              
      else 
!-- Cooperstein, Van den Horn, Baron method                             
         prate=fp4 
      end if 
      prate=c2cu*prate/beta5 
!                                                                       
!--if nu trapped, bail out without computing energy loss                
      if(iflag.eq.1) then 
         due=0. 
         dup=0. 
         return 
      endif 
!                                                                       
!--compute cooling due to electron capture (per nucleon)                
!                                                                       
      if (icvb.eq.1) then 
         due=fe5 + p1*fe4*beta + p2*fe3*beta2 + p3*fe2*beta3 +          &
     &        p4*fe1*beta4 + p5*fe0*beta5                               
      else 
!-- Cooperstein, Van den Horn, Baron method                             
         due=fe5 
      end if 
      due=c3cu*due/beta6 
!                                                                       
!--compute cooling due to positron capture (per nucleon)                
!                                                                       
      if (icvb.eq.1) then 
         dup=fp5 - p1*fp4*beta + p2*fp3*beta2 - p3*fp2*beta3 +          &
     &        p4*fp1*beta4 - p5*fp0*beta5                               
      else 
         dup=fp5 
      end if 
      dup=c3cu*dup/beta6 
!                                                                       
      return 
      END                                           
!                                                                       
      subroutine step(ntstep) 
!************************************************************           
!                                                           *           
!  This subroutine integrates the system of equations using *           
!  a Runge-Kutta-Fehlberg integrator of second order.       *           
!  Particles are allowed to have individual time-steps.     *           
!  All particles are synchronized every dtime at which time *           
!  the subroutine is exited.                                *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx, ntstep, ind, mlin_grid_size
      logical from_dump, post_bounce
      character*128 mlmodel_name, rho_file, x_file      
      character*10 frmtx
      character*11 frmtrho
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
!                                                                       
      logical ifirst, reset(idim)
!                                                                       
      logical ebetaeq, pbetaeq 
      common /beta/ ebetaeq(idim), pbetaeq(idim) 
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /carac/ deltam(idim), abar(idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /cgas / gamma 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
      common /core / dcore, xmcore 
      common /ener2/ tkin, tterm 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /epcap/ betafac, c2cu, c3cu 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /freez/ ufreez(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /numb / ncell, ncell1 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)                  
      common /prev / xold(0:idim),vold(0:idim),rhold(idim),             &
     &              prold(idim), tempold(idim), yeold(idim),            &
     &              xnold(idim), xpold(idim)                            
      common /shock/ cq,cl 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /swesty/ xpf(idim), pvar2(idim), pvar3(idim), pvar4(idim) 
      common /tempe/ temp(idim) 
      common /therm/ xmu(idim) 
      common /timej / time, dt
      common /rshock/ shock_ind, shock_x
      common /pns/ pns_ind, pns_x

      logical trapnue, trapnueb, trapnux, print_nuloss 
      common /trap / trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /typef/ iextf, ieos 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /dum  / dumx(0:idim),dumv(0:idim),dumu(idim),              &
     &               dumye(idim)                                        
      common /dum2 / dumynue(idim),dumynueb(idim),dumynux(idim),        &
     &               dumunue(idim),dumunueb(idim),dumunux(idim)         
      common /f1   / f1v(0:idim),f1u(idim),f1ye(idim) 
      common /f2   / f2v(0:idim),f2u(idim),f2ye(idim) 
      !common /fturb/ geff(idim), fmix1(idim), fmix2(idim) 
      !common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
      !common /dturb/ dumvt2(idim) 
      common /f1nu / f1ynue(idim),f1ynueb(idim),f1ynux(idim),           &
     &               f1unue(idim),f1unueb(idim),f1unux(idim)            
      common /f2nu / f2ynue(idim),f2ynueb(idim),f2ynux(idim),           &
     &               f2unue(idim),f2unueb(idim),f2unux(idim)            
      common /nuout/ rlumnue, rlumnueb, rlumnux,                        &
     &               enue, enueb, enux, e2nue, e2nueb, e2nux            
!                                                                       
      common /timei/ istep(idim),t0(idim),steps(idim),                  &
     &               dum2v(idim)                                        
      common /propt/ dtime,tmax 
      common /nsat/ satc,xtime 
      common /nuprint/ nups,nupk,tacr 
      common /outp/ rout, p1out, p2out
      common /mlmod/ mlmodel_name
      common /dump/ from_dump
      common /bnc/ rlumnue_max, bounce_ntstep, bounce_time, post_bounce
      common /interp/ mlin_grid_size

      save ifirst
      data ifirst/.true./      
      data tiny/3.e-5/         
!                                                                       
!--Compute next dump time and initialise variables.                     
!                                                                       
      tol=1.d-3 
      dt0=1.d0 
      tnext=time+dtime 
!      xlog2=0.30103                                                    
!                                                                       
!--define coefficients for Runge-Kutta integrator                       
!                                                                       
      f11=0.5*256./255. 
      f21=1./256. 
      f22=255./256. 
      e1=1./512. 
!                                                                       
!--first time set time step to default                                  
!
      !print*, steps
      !if ( ANY( steps>(dt0/1.d8) ) ) then
      !    from_dump=.true.
      !    print*, '-----------> From Dump is True!', steps(1)
      !endif            

      if(ifirst) then
         ifirst=.false. 
         dumx(0)=x(0) 
         do i=1,ncell
            if (from_dump.eqv..false.) then
                steps(i)=dt0/1.d8
            end if
            istep(i)=nint(dtime/steps(i)) 
            t0(i)=time 
            tempold(i)=temp(i) 
            rhold(i)=rho(i) 
            yeold(i)=ye(i) 
            xnold(i)=xn(i) 
            xpold(i)=xp(i) 
         enddo 
         
         write(*,'(A)')'<      Hydro Setup     >'       
         print_nuloss=.false. 
         call hydro(time,ncell,x,v,                                     &
     &           u,rho,ye,f1v,f1u,f1ye,q,                               &
     !fmix1,                         &
     &           ynue,ynueb,ynux,f1ynue,f1ynueb,f1ynux,                 &
     &           unue,unueb,unux,f1unue,f1unueb,f1unux,                 &
     &           print_nuloss,ntstep)                                          
      end if 
      !do i=30,40
      !    print*, i, ye(i)
      !enddo      
!                                                                       
!--set step counter                                                     
!                                                                       
      do i=1,ncell 
         istep(i)=nint(2.*dtime/steps(i)) 
      enddo 
!                                                                       
!--integration                                                          
!  -----------                                                          
!                                                                       
!  a) Unique time step case                                             
!  ------------------------                                             
!                                                                       
      ivar=0 
      if(ivar.eq.0) then 
!                                                                       
!--get predictions at half time step                                    
!--start of the new time step - essentially this is the main loop       
!                                                                       
   99    continue 
!                                                                       
!--gas particles                                                        
!                                                                       
         !print *, '[steps] ',steps(1),x(0),x(1)                        
         do i=1,ncell 
            dtf11=f11*steps(i) 
            dumx(i)=x(i) + dtf11*v(i) 
            dumv(i)=v(i) + dtf11*f1v(i) 
            if (time.lt.0.02.and.rho(i).lt.5.d6)                        &
     &           dumv(i)=min(dumv(i),-0.1)                              
            dumu(i)=u(i) + dtf11*f1u(i) 
            dumye(i)=ye(i) + dtf11*f1ye(i) 
            dumynue(i)=ynue(i) + dtf11*f1ynue(i) 
            dumynueb(i)=ynueb(i) + dtf11*f1ynueb(i) 
            dumynux(i)=ynux(i) + dtf11*f1ynux(i) 
            dumunue(i)=unue(i) + dtf11*f1unue(i) 
            dumunueb(i)=unueb(i) + dtf11*f1unueb(i) 
            dumunux(i)=unux(i) + dtf11*f1unux(i) 
            !dumvt2(i)=vturb2(i) + dtf11*fmix1(i) 
            reset(i)=.false. 
            if (dumye(i).lt.0.02) then 
               dumye(i)=.02 
            end if 
         enddo 
!                                                                       
!--get forces at half the time step, using predictions                  
!                                                                       
         thalf=time + f11*steps(1) 
         print_nuloss=.false. 
         !call write_data(ntstep*2,ncell,x,v)

         call hydro(thalf,ncell,dumx,dumv,                              &
     &         dumu,rho,dumye,f2v,f2u,f2ye,q,                           &
     !fmix2,                     &
     &         dumynue,dumynueb,dumynux,f2ynue,f2ynueb,f2ynux,          &
     &         dumunue,dumunueb,dumunux,f2unue,f2unueb,f2unux,          &
     &         print_nuloss,ntstep)        
!                                                                       
!--advance all the gas particles                                        
!        
         !do i=30,40
         !    print*, i, dtf21*f1ye(i), dtf22*f2ye(i) 
          !print*, i, dumye(i), ye(i),dtf21*f1ye(i), dtf22*f2ye(i) 
          !print*, i, ye(i)+dtf21*f1ye(i)+dtf22*f2ye(i) 
         !enddo
         do i=1,ncell 
            dtf21=f21*steps(i) 
            dtf22=f22*steps(i) 
            dumx(i)=x(i) + dtf21*v(i) + dtf22*dumv(i) 
            dumv(i)=v(i) + dtf21*f1v(i) + dtf22*f2v(i) 
            if (time.lt.0.02.and.rho(i).lt.5.d6)                        &
     &           dumv(i)=min(dumv(i),-0.1)                              
            dumu(i)=u(i) + dtf21*f1u(i) + dtf22*f2u(i) 
            dumye(i)=ye(i) + dtf21*f1ye(i) + dtf22*f2ye(i) 
            dumynue(i)=ynue(i) + dtf21*f1ynue(i) + dtf22*f2ynue(i) 
            dumynueb(i)=ynueb(i) + dtf21*f1ynueb(i) + dtf22*f2ynueb(i) 
            dumynux(i)=ynux(i) + dtf21*f1ynux(i) + dtf22*f2ynux(i) 
            dumunue(i)=unue(i) + dtf21*f1unue(i) + dtf22*f2unue(i) 
            dumunueb(i)=unueb(i) + dtf21*f1unueb(i) + dtf22*f2unueb(i) 
            dumunux(i)=unux(i) + dtf21*f1unux(i) + dtf22*f2unux(i) 
            !dumvt2(i)=vturb2(i) + dtf21*fmix1(i) + dtf22*fmix2(i) 
            if (dumye(i).lt.0.02) then 
               dumye(i)=.02 
            end if 
         enddo 
         !do i=30,40
         !    print*, i, dft21, dft22, steps(i)
         !enddo         
!                                                                       
!     set saturation const=0                                            
         satc=0 
!                                                                       
!--get forces at end of time step                                       
!                                                                       
         tfull=time + steps(1) 
         print_nuloss=.true. 
         !call write_data(ntstep*2+1,ncell,x,v)
         !print*, '.......................'
         !do i=30,40
         ! print*, i, dumye(i)
         !enddo         
         !print*, '-----------------------'
         call hydro(tfull,ncell,dumx,dumv,                              &
     &         dumu,rho,dumye,f2v,f2u,f2ye,q,                           &
     !fmix2,                     &
     &         dumynue,dumynueb,dumynux,f2ynue,f2ynueb,f2ynux,          &
     &         dumunue,dumunueb,dumunux,f2unue,f2unueb,f2unux,          &
     &         print_nuloss,ntstep)    
      !call exit(0)     
!                                                                       
!--estimate integration error and set time step. If reduction           
!  the maximum time step is set to dtime.                               
!  Error checking is supressed for particles having been                
!  reset.                                                               
!                                                                       
         rapmax=0. 
         ermax=-1.e10 
         stepmin=1.e30 
         stepmax=0. 
         ierr=0 
         nreset=0 
         denmax=0. 
         tempmax=0. 
         do i=1,ncell 
            if(reset(i)) then 
               steps(i)=1.e20 
               nreset=nreset + 1 
            else 
               dxi=(x(i)-x(i-1)) 
               vsi=0.05*vsound(i) 
               ysi=ye(i) 
               if (ifleos(i).eq.1) then 
                  usi=u(i) 
               elseif (ifleos(i).eq.2) then 
!--this because in NSE, we can have u(i)< lt 0                          
!--so take energy with respect to iron                                  
                  usi=u(i)+8.8*umevnuc 
               elseif (ifleos(i).eq.3.or.ifleos(i).eq.4) then
                  if (ieos.eq.3) then 
                     usi=u(i)+8.8*umevnuc 
                  else 
!--variable of state is entropy                                         
                     usi=u(i)*temp(i) 
                  endif 
               else
                  usi=u(i)                   
               endif 
               erx=abs(dumv(i))/dxi 
               erv=abs(f1v(i)-f2v(i))/(abs(v(i))+vsi) 
               eru=abs(f1u(i)-f2u(i))/u(i) 
               erye=abs(f1ye(i)-f2ye(i))/ysi 
               erynue=abs(f1ynue(i)-f2ynue(i))/(ysi+ynue(i)) 
               erynueb=abs(f1ynueb(i)-f2ynueb(i))/(ysi+ynueb(i)) 
               erynux=abs(f1ynux(i)-f2ynux(i))/(ysi+ynux(i)) 
               erunue=abs(f1unue(i)-f2unue(i))/(usi+unue(i)) 
               erunueb=abs(f1unueb(i)-f2unueb(i))/(usi+unueb(i)) 
               erunux=abs(f1unux(i)-f2unux(i))/(usi+unux(i)) 
!                                                                       
!--get maximum error and determine time step                            
!                                                                       
               erm=dmax1(erx,erv,eru,erye,                              &
     &            erynue,erynueb,erynux,erunue,erunueb,erunux)          
               if(erm.gt.ermax)then 
                  if(erm.eq.erx)then 
                     ierr=1 
                  elseif(erm.eq.erv)then 
                     ierr=4 
                  elseif(erm.eq.eru)then 
                     ierr=7 
                  elseif(erm.eq.erye) then 
                     ierr=9 
                  elseif(erm.eq.erynue) then 
                     ierr=10 
                  elseif(erm.eq.erynueb) then 
                     ierr=11 
                  elseif(erm.eq.erynux) then 
                     ierr=12 
                  elseif(erm.eq.erunue) then 
                     ierr=13 
                  elseif(erm.eq.erunueb) then 
                     ierr=14 
                  elseif(erm.eq.erunux) then 
                     ierr=15 
                  endif 
                  imax=i 
                  ermax=erm 
               endif 
               rap=steps(i)*erm*e1/tol + tiny 
               rapmax=dmax1(rapmax,rap) 
               rmod=min(2.0d0,1./dsqrt(rap)) 
               steps(i)=min(steps(i)*rmod,dtime,dt0) 
!               if (time.gt.1.d-3) then                                 
!                  print *, steps(i)*rmod,steps(i),rmod,dtime,dt0       
!               end if                                                  
               tempmax=max(temp(i),tempmax) 
               denmax=max(rho(i),denmax) 
            end if 
         enddo 
!         write(*,*)'max error flag: ',ierr,imax,ebetaeq(imax)          
!         write(*,*)'temp,x',temp(imax),dumx(imax),dumx(imax-1)         
!         write(*,*)'max temperature, max density', tempmax,denmax      
!         write(*,*)'v,rho',dumv(imax),rho(imax)                        
!         write(*,*)'u,ye',u(imax),ye(imax)                             
!         write(*,*)'new: u,ye',dumu(imax),dumye(imax)                  
!         write(*,*)'fx',f2v(imax)                                      
!         write(*,*)'du,dye',f2u(imax),f2ye(imax)                       
!         write(*,*)'ynue,ynueb,ynux',                                  
!     1              dumynue(imax),dumynueb(imax),dumynux(imax)         
!         write(*,*)'dynue,dynueb,dynux',                               
!     1              f2ynue(imax),f2ynueb(imax),f2ynux(imax)            
!         write(*,*)'unue,unueb,unux',                                  
!     1              dumunue(imax),dumunueb(imax),dumunux(imax)         
!         write(*,*)'dunue,dunueb,dunux',                               
!     1              f2unue(imax),f2unueb(imax),f2unux(imax)            
!                                                                       
!--find minimum time step                                               
!                                                                       
         dtmin=1.e30 
         do i=1,ncell 
            dtmin=dmin1(dtmin,steps(i)) 
            stepmax=dmax1(stepmax,steps(i)) 
            stepmin=dmin1(stepmin,steps(i)) 
         enddo 
         xtime=xtime*0.8d0**satc 
         !print *,'xtime=',xtime,'satc=',satc                           
         do i=1,ncell 
            steps(i)=min(dtmin,xtime) 
            if (time.lt.1.d-7) then 
               steps(i)=min(steps(i),1.d-9) 
            end if 
         enddo 
!         write(*,110)'step:t,dt,rapmax,nreset', tfull,dtmin,rapmax,    
!     1             nreset                                              
  110    format(A30,3(1g12.5,1x),I3) 
!                                                                       
!--check for rejection                                                  
!                                                                       
         !print *, 'XXXXXXXXXX rapmax',rapmax                           
         if(rapmax.gt.1.5)then 
            go to 99 
         end if 
!                                                                       
!--update system                                                        
!                                                                       
         !print *, 'XXXXXXXXXXXXXXXXXXXX update system'                 
         do i=1,ncell 
            x(i)=dumx(i) 
            v(i)=dumv(i) 
            u(i)=dumu(i) 
            ye(i)=dumye(i) 
            f1v(i)=f2v(i) 
            f1u(i)=f2u(i) 
            f1ye(i)=f2ye(i) 
            ynue(i)=dumynue(i) 
            ynueb(i)=dumynueb(i) 
            ynux(i)=dumynux(i) 
            f1ynue(i)=f2ynue(i) 
            f1ynueb(i)=f2ynueb(i) 
            f1ynux(i)=f2ynux(i) 
            unue(i)=dumunue(i) 
            unueb(i)=dumunueb(i) 
            unux(i)=dumunux(i) 
            !vturb2(i)=dumvt2(i) 
            f1unue(i)=f2unue(i) 
            f1unueb(i)=f2unueb(i) 
            f1unux(i)=f2unux(i) 
         enddo 
         if (x(ncell).gt.rout) then 
            pr(ncell+1)=2.*p1out 
         else 
            pr(ncell+1)=0.1*p1out 
         end if 
!                                                                       
!--burning                                                              
!                                                                       
         iburn=0 
         !print *, 'xxxxxx',iburn                                       
         if(iburn.eq.1)then 
            iflg=0 
            !print *, 'call burn - time',time,steps(1)                  
            call burn(ncell,iflg,tfull,steps(1)) 
            if(iflg.eq.1)then 
               call eosflg(ncell,rho,ye,u,f1ye,f1u) 
               do i=1,ncell 
                  tempold(i)=temp(i) 
                  rhold(i)=rho(i) 
                  yeold(i)=ye(i) 
                  xnold(i)=xn(i) 
                  xpold(i)=xp(i) 
               end do 
               call hydro(tfull,ncell,x,v,                              &
     &              u,rho,ye,f1v,f1u,f1ye,q,                            &
     !fmix1,                      &
     &              ynue,ynueb,ynux,f1ynue,f1ynueb,f1ynux,              &
     &              unue,unueb,unux,f1unue,f1unueb,f1unux,              &
     &              print_nuloss,ntstep)                                       
            endif 
         endif 
         time=tfull 
!                                                                       
!--remove particles if necessary                                        
!                                                                       
!         call rmp(tfull,ncell,x,v,u,rho,ye,rb,                         
!     1            ynue,ynueb,ynux,unue,unueb,unux,nrem)                
!         if (time.gt.tacr) then                                        
!            tacr=tacr+2.d-4                                            
!            call amp(tfull,ncell,x,v,u,rho,ye,rb,                      
!     1           ynue,ynueb,ynux,unue,unueb,unux)                      
!         end if                                                        
!         print *,'rmp',xmcore,ncell,rb,tacr                            
!         if (nrem.gt.0)then                                            
!            call hydro(tfull,ncell,x,v,                                
!     $           u,rho,ye,f1v,f1u,f1ye,q,fmix1,                        
!     $           ynue,ynueb,ynux,f1ynue,f1ynueb,f1ynux,                
!     $           unue,unueb,unux,f1unue,f1unueb,f1unux)                
!            do i=1,ncell                                               
!               steps(i)=steps(i)*1.d-1                                 
!            end do                                                     
!         end if                                                        
!                                                                       
!--flag particles according to physical state                           
!
         if (ieos.ne.5) then                                                                       
            call eosflg(ncell,rho,ye,u,f1ye,f1u) 
         endif
!                                                                       
!--do various neutrino flagging and checking                            
!--skip neutrino                                                        
         call nucheck(ncell,x,ye,                                       &
     &            ynue,ynueb,ynux,unue,unueb,unux)                      
!                                                                       
         do i=1,ncell 
            tempold(i)=temp(i) 
            rhold(i)=rho(i) 
            yeold(i)=ye(i) 
            xnold(i)=xn(i) 
            xpold(i)=xp(i) 
         enddo 
         convf=ufoe/utime 
         nupk=nupk+1 
!         if (time.gt.4.d-3) then                                       
!            5000                                                       
!         end if                                                        
         if (mod(nupk,nups).eq.0) then 
            write(59,120)time,convf*rlumnue,enue,convf*rlumnueb,enueb,  &
     &           convf*rlumnux,enux                                     
            ke=0. 
            do i=1,ncell 
               if (v(i).gt.0) then 
                  ke=ke+deltam(i)*v(i)**2 
               end if 
            end do 
!            write(58,120) time*10., ke*2.d-2                           
!            write(90,120)time,x(1),x(2),x(3),x(4),x(5),x(6)            
!            write(91,120)time,x(7),x(8),x(9),x(10),x(11),x(12)         
!            write(92,120)time,x(13),x(14),x(15),x(16),x(17),x(18)      
!            write(93,120)time,x(19),x(20),x(21),x(22),x(23),x(24)      
!            write(94,120)time,x(25),x(26),x(27),x(28),x(29),x(30)      
!            write(95,120)time,x(32),x(34),x(36),x(38),x(40),x(42)      
!            write(96,120)time,x(44),x(46),x(48),x(50),x(55),x(60)      
!            write(97,120)time,x(65),x(70),x(75),x(80),x(85),x(90)      
!            write(98,120)time,x(100),x(110),x(120),x(130),x(140),x(150)
!            write(99,120)time,x(160),x(170)                            
         end if 
!                                                                       
  120    format((1pe12.5),6(1x,1pe10.2))
!                                                                       
!--start new time step                                                  
!                                                                       
         if(time.lt.tnext)then 
            !print *, '****time*****',time                              
  520       format(A,I12,A) 
            print 520,'<',ntstep,                                       &
     &      '          > ----------------------------------'            
            write(*,500)'[    time/tmax, dt (s) ]',                     &
     &           time*10,'/',tmax*10,steps(1)*10                              
  500       format(A,1p,E10.3,A,E9.3,E15.3)
  
            if (post_bounce.eqv..true.) then
                write(*,501)'[    bounce time (s)   ]',                     &
         &           bounce_time*10                              
      501       format(A,1p,E10.3)
            endif
!
            if (mod(nupk,nups).eq.0) then 
               lu=72 
               open(lu,file='1dtmp',form='unformatted') 
               call printout(lu) 
               close (lu) 
            end if 
            ntstep=ntstep+1 
            go to 99 
         end if 
         return 
      end if 
!                                                                       
      END  


      subroutine write_data(ntstep, ncell,x, v)
       
      implicit double precision (a-h,o-z) 

      dimension x(0:ncell),v(0:ncell) 
      !integer :: i,j,shock_ind, ind
      character*30 :: rho_file, x_file
      
      
      if (ntstep<10) then          
          write (rho_file, "(A9,I1,A4)") "tmpdir/v_", ntstep, ".txt"
          write (x_file, "(A9,I1,A4)") "tmpdir/x_", ntstep, ".txt"         
      else if (ntstep<100) then
          write (rho_file, "(A9,I2,A4)") "tmpdir/v_", ntstep, ".txt"
          write (x_file, "(A9,I2,A4)") "tmpdir/x_", ntstep, ".txt"
      else if (ntstep<1000) then
          write (rho_file, "(A9,I3,A4)") "tmpdir/v_", ntstep, ".txt"
          write (x_file, "(A9,I3,A4)") "tmpdir/x_", ntstep, ".txt"
      else if (ntstep<10000) then
          write (rho_file, "(A9,I4,A4)") "tmpdir/v_", ntstep, ".txt"
          write (x_file, "(A9,I4,A4)") "tmpdir/x_", ntstep, ".txt"
      else
          write (rho_file, "(A9,I5,A4)") "tmpdir/v_", ntstep, ".txt"
          write (x_file, "(A9,I5,A4)") "tmpdir/x_", ntstep, ".txt"
      end if 

      if (MOD(ntstep, 1).eq.0) then
          open(unit=995, file=trim(rho_file), status='replace')
          do ind=1,ncell
            write(995,*), v(ind)
          end do
          close(unit=995)

          open(unit=996, file=trim(x_file), status='replace')
          do ind=1,ncell
            write(996,*), x(ind)
          end do
          close(unit=996)
          !print*, '=== end rho ==='
          !call EXIT(0)
      end if
      end
      
      subroutine unit 
!************************************************************           
!                                                           *           
!  this routine computes the transformation between the     *           
!  physical units (cgs) and the units used in the code.     *           
!  And the value of physical constants in code units        *           
!                                                           *           
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      double precision umass 
      double precision usltemp, uslrho, uslu, uslp,u2slu 
      double precision ud1,ut1,ue1,ueg1,uec1 
      double precision gg1, arad1, bigr1 
      double precision uopr, uotemp, uorho1, uotemp1, uou1 
      double precision ufoe1, sigma1, sigma2, xsecnn1, xsecne1,fermi 
!                                                                       
      common /typef/ iextf, ieos 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu 
      common /epcap/ betafac, c2cu, c3cu 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne 
!                                                                       
      data ggcgs/6.67e-8/, avo/6.02e23/ 
      data aradcgs /7.565e-15/, boltzk/1.381e-16/ 
      data hbar/1.055e-27/, ccgs/3e10/ 
      data emssmev/0.511e0/, boltzmev/8.617e-11/ 
      data ergmev /6.2422e5/, sigma1/9d-44/, sigma2/5.6d-45/ 
      data c2cgs /6.15e-4/, c3cgs /5.04e-10/, fermi/1d-13/ 
!                                                                       
!--1) work out code units:                                              
!                                                                       
!--specifie mass unit (g)                                               
!                                                                       
      umass=2d33 
!                                                                       
!--specifie distance unit (cm)                                          
!                                                                       
      udist=1e9 
!                                                                       
!--transformation factor for :                                          
!                                                                       
! 1a) density                                                           
!                                                                       
      ud1=dble(umass)/dble(udist)**3 
      udens=ud1 
!                                                                       
! 1b) time                                                              
!                                                                       
!     ut1=dsqrt(dble(udist)**3/(dble(gg)*umass))                        
!     utime=ut1                                                         
      ut1=1.d1 
      utime=ut1 
!                                                                       
! 1c) ergs                                                              
!                                                                       
      ue1=dble(umass)*dble(udist)**2/dble(utime)**2 
!--uerg will overflow                                                   
!      uerg=ue1                                                         
!                                                                       
! 1d) ergs per gram                                                     
!                                                                       
      ueg1=dble(udist)**2/dble(utime)**2 
      uergg=ueg1 
!                                                                       
! 1e) ergs per cc                                                       
!                                                                       
      uec1=dble(umass)/(dble(udist)*dble(utime)**2) 
      uergcc=uec1 
!                                                                       
!--2) constants                                                         
!                                                                       
! 2a) gravitational                                                     
!                                                                       
      gg1=dble(ggcgs)*umass*dble(utime)**2/(dble(udist)**3) 
      gg=gg1 
!                                                                       
! 2aa) velocity of light                                                
!                                                                       
      clight=dble(ccgs*utime/udist) 
                                                                        
!                                                                       
! 2b) Stefan Boltzmann (note that T unit is 1e9 K here)                 
!                                                                       
      utemp=1e9 
      arad1=dble(aradcgs)/ue1*dble(udist)**3*dble(utemp)**4 
      arad=arad1 
!                                                                       
! 2c) Perfect Gas "R" constant                                          
!                                                                       
      bigr1=dble(avo)*dble(boltzk)/ue1*umass*dble(utemp) 
      bigr=bigr1 
!                                                                       
! 2d) nucleon+nu x-section/4pi, in udist**2/umass/Mev**2                
!                                                                       
      xsecnn1=sigma1*umass*dble(avo)/dble(udist*udist) 
      xsecnn=xsecnn1/(4d0*3.14159d0) 
!                                                                       
! 2e) e+nu x-section/4pi, in Enu(Mev)*udist**2/umass/Mev**2/utemp**4    
!                                                                       
      xsecne1=sigma2*umass*dble(ergmev)*(dble(utemp)*dble(boltzk))**4/  &
     &     (dble(udist)**2*ud1*3.14159d0*(dble(hbar)*dble(ccgs))**3)    
      xsecne=xsecne1/(4d0*3.14159d0) 
!                                                                       
!--3a) Conversion to the Ocean eos units from code units:               
!     in the ocean eos; density unit=1e7g/cc                            
!                       temperature unit= 1e9K                          
!                       energy/mass=1.0e17cgs                           
!                       energy/vol =1.0e24cgs                           
!                                                                       
      uotemp=1d9/dble(utemp) 
      uotemp1=1.d0/uotemp 
      uopr=1.d24/dble(uergcc) 
      uorho1=ud1/1.d7 
      uou1=dble(uergg)/1.d17 
                                                                        
                                                                        
                                                                        
                                                                        
!                                                                       
!--3b) Conversion to the Swesty-Lattimer units from code units:         
!                                                                       
      usltemp=utemp*boltzmev 
      uslrho=ud1*avo*fermi**3 
      uslp=uec1*fermi**3/1.602d-6 
      if (ieos.eq.3) then 
!--if u is internal energy                                              
         uslu=dble(ergmev)*dble(uergg)/dble(avo) 
         u2slu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk)) 
      elseif (ieos.eq.4) then 
!--if u is specific entropy                                             
         uslu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk)) 
         u2slu=dble(ergmev)*dble(uergg)/dble(avo) 
      endif 
      print *,'ieos,uslu',ieos,uslu 
!                                                                       
!--4) common unit2 stuff                                                
!                                                                       
! 4a) temp code unit in Mev                                             
!                                                                       
      utmev=utemp*boltzmev 
!                                                                       
! 4b) energy code unit in foes                                          
!                                                                       
      ufoe1=ue1/1d51 
      ufoe=ufoe1 
!                                                                       
! 4c) Mev/nucleon in code units                                         
!                                                                       
      umevnuc=avo/(ergmev*uergg) 
!                                                                       
! 4e) Mev to errgs times avogadro's number                              
!                                                                       
      umeverg=avo*1.602e-6 
!                                                                       
!--5) epcapture betafac, c2, and c3.                                    
!                                                                       
      betafac=emssmev/utmev 
      c2cu=c2cgs*utime 
      c3cu=c3cgs*utime*ergmev 
!                                                                       
      write(*,100)umass,udist,udens,utime,uergg,uergcc 
  100 format(1x,5(1pe12.5,1x),/,1x,2(1pe12.5,1x)) 
!                                                                       
   99 return 
      END                                           
                                                                        
      subroutine burn(ncell,iflg,time,deltat) 
!*****************************************************************      
!                                                                *      
!  Nuclear network subroutine.                                   *      
!  This subroutine uses a 14 elements alpha network from fkt     *      
!                                                                *      
!*****************************************************************      
!                                                                       
      parameter (idim=10000) 
      parameter (idim1=idim+1) 
      parameter (iqn=17) 
      parameter (nel=14) 
!                                                                       
      integer ifleos 
      double precision pf(7,17),pb(7,14),rrat(17),rlam(15),ya(15),      &
     &                 def(14),ff(14),yold(15),t9,rhoi,tburn,dtold,     &
     &                 dtburn,dfmin,dfm,err,eps,dtsec,enrg,             &
     &                 errmax,scale,a(14),yinput(15),zn(14),            &
     &                 yel,coefg,coeft,ar(14),ar11,ui,tempi,            &
     &                 abari,yei,rhoit,xpn,etai,time                    
                                                                        
      double precision deltam,abar,u,rho,ye,q,xp,xn,eta,xalpha,         &
     &     xheavy,yeh,pr,vsound,u2,vsmax,temp,amas,znum                 
      double precision umass,udist,udens,utime,                         &
     &     uergg,uergcc,gamma,uvit,uit,deltat                           
      double precision t9nse, rhoswe, rhonue, rhonux, xmue,             &
     &     xmuhat                                                       
      double precision ptot,cs,yehi,stot,ufact,yefact 
      double precision ufreez 
!                                                                       
      common /ceos / amas(iqn), znum(iqn) 
      common /cc   / ycc(idim,iqn), yccave(iqn) 
      common /carac/ deltam(idim), abar(idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /hnucl/ xalpha(idim),xheavy(idim), yeh(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      logical ifign 
      common /ign  / ifign(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /cgas/ gamma 
      common /cases/ t9nse, rhoswe, rhonue, rhonux 
      common /cpots/ xmue(idim), xmuhat(idim) 
      common /freez/ ufreez(idim) 
      common /ufactor/ ufact,yefact 
                                                                        
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3 
      common /nucen/ totnuc, stepnuc 
      common /tempe/ temp(idim) 
                                                                        
      common /rateinfo/ ratforward(7,17),ratbackward(7,14) 
!                                                                       
      logical first 
!                                                                       
      character*2 element(nel) 
!                                                                       
      data   amas /  1.d00,  1.0d0, 1.0d0, 4.d00, 12.d00, 16.d00,       &
     &              20.d00, 24.d00,28.d00, 32.d00, 36.d00, 40.d00,      &
     &              44.0d0, 48.0d0, 52.0d0, 56.0d0, 60.0d0/             
      data   znum /  0.d00,  1.0d0,  0.0d0, 2.d00,  6.d00,  8.d00,      &
     &              10.d00, 12.d00, 14.d00, 16.d00, 18.d00, 20.d00,     &
     &              22.0d0, 24.0d0, 26.0d0, 28.0d0, 30.0d0/             
!                                                                       
      data eps/1.d-3/ , maxstep/10000/ 
      data first/.true./ 
      data a/4.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,            &
     &       56.,60./                                                   
      data zn/2.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30./ 
      data element/'He','C ','O ','Ne','Mg','Si','S ','Ar','Ca','Ti',   &
     &             'Cr','Fe','Ni','Zn'/                                 
      nqn=17 
      idisk2=10 
!                                                                       
!--if first passage read rates                                          
!                                                                       
      if(first) then 
!         write(*,100)nel,(i,element(i),i=1,nel)                        
!  100    format(/,' b u r n   m o d e l  : (F. K. Thielemann)',/,      
!     1            ' -------------------',//,                           
!     1            ' Total number of elements considered : ',i4,/,      
!     2   4(1x,i2,' : ',a2,3x,i2,' : ',a2,3x,i2,' : ',a2,3x,i2,' : ',   
!     3   a2,/))                                                        
!         write(*,101)                                                  
!  101    format(/,'  reactions considered : ',/)                       
         open(idisk2,file='alphanet.dat') 
         call rrate(pf,pb,iprint,idisk2) 
         do i=1,7 
            do j=1,17 
               ratforward(i,j)=pf(i,j) 
            end do 
         end do 
         do i=1,7 
            do j=1,14 
               ratbackward(i,j)=pb(i,j) 
            end do 
         end do 
                                                                        
!         write(*,102)                                                  
!  102    format(///)                                                   
         first=.false. 
         close(idisk2) 
         open(80,file='energy.nuc') 
      end if 
!                                                                       
         do i=1,7 
            do j=1,17 
               pf(i,j)=ratforward(i,j) 
            end do 
         end do 
         do i=1,7 
            do j=1,14 
               pb(i,j)=ratbackward(i,j) 
            end do 
         end do 
                                                                        
!                                                                       
!--initialisation                                                       
!                                                                       
      stepnuc=0.0 
      dtsec=deltat*utime 
      uvit=udist/utime 
      errmax=0. 
      input=2 
!                                                                       
!--compute the burn for all particles                                   
!                                                                       
      do icell=1,ncell 
         if(.not.ifign(icell))go to 95 
         t9=temp(icell) 
         if(t9.lt.0.08)go to 95 
!         if(t9.lt.0.12)go to 95                                        
         iflg=1 
         rhoi=rho(icell)*udens 
         ui=u(icell)*uergg 
         do iel=1,nel 
            yold(iel)=0.d0 
            ya(iel)=0.d0 
            ya(iel)=dble(ycc(icell,iel+3)) 
            yinput(iel)=ya(iel) 
         enddo 
         tburn=0.0 
!                                                                       
!--begin of new time step ka                                            
!                                                                       
         do ka=1,maxstep 
            dtold=dtburn 
!                                                                       
!--calculate rates for temperature                                      
!                                                                       
            call genpar(rhoi,t9,a,zn,ya,yel,ar,ar11,coefg,coeft) 
            call rates(icell,t9,pf,pb,rrat,rlam,yel,coefg,coeft,a,ar,   &
     &                 ar11,zn)                                         
!                                                                       
!--compute ydot(i)=ff(i) : system of differential equations             
!  with present ya(i) rho and rates                                     
!                                                                       
            call derivn(ya,rhoi,ff,rrat,rlam) 
!                                                                       
!--determine time step from dt=0.1*min (ya(i)/ydot(i))                  
!  for first time step (ka=1) take ydot=ff                              
!  for further time steps take ydot=(ya-yold)/dtold                     
!                                                                       
            dfmin=1.d20 
            if(ka.ne.1) goto 60 
            do 50 i=1,nel 
               if(ya(i).lt.eps.or.ff(i).eq.0.0d0) goto 50 
               dfm=abs(ya(i)/ff(i)) 
               dfmin=min(dfm,dfmin) 
   50       continue 
            goto 80 
   60       do 70 i=1,nel 
               def(i)=(ya(i)-yold(i)) 
               if(ya(i).lt.eps.or.def(i).eq.0.0) goto 70 
               dfm=abs(ya(i)/def(i))*dtold 
               dfmin=min(dfmin,dfm) 
   70       continue 
   80       dtburn=7.0d-2*dfmin 
            dtbold=dtburn 
            dtburn=min(dtburn,dtsec-tburn) 
!                                                                       
!--reset yold(i) by updating from last time step                        
!                                                                       
   85       do i=1,nel 
               yold(i)=ya(i) 
            enddo 
!                                                                       
!--calculate new abundances at t+dt                                     
!                                                                       
   35       continue 
            call newab(yold,ya,rhoi,rrat,rlam,ff,dtburn,err,eps) 
            if(abs(err).gt.abs(errmax))errmax=err 
            tburn=tburn+dtburn 
            call epsb(yold,ya,nel,enrg) 
!                                                                       
!--if energy has been generated, compute new T                          
!                                                                       
            ui = ui + enrg 
            uit=ui/uergg 
            rhoit=rhoi/udens 
            tempi=t9 
            yei=yel 
            abari=0.d0 
            do j=4,nqn 
               abari=abari+dble(ya(j-3))*amas(j)**2 
            enddo 
            if(tburn.ge.dtsec)then 
               iout=1 
               call rootemp2(icell,rhoit,uit,tempi,yei,abari,           &
     &              ptot,cs,etai,stot)                                  
               t9=tempi 
               pr(icell)=ptot/uergcc 
               vsound(icell)=cs/uvit 
               abar(icell)=abari 
               go to 92 
            else 
               iout=0 
               if (tempi.gt.t9nse+1.) then 
                  ufreez(icell)=-ufact*8.7d2 
                  temp(icell)=tempi 
                  goto 95 
               else 
                  call rootemp2(icell,rhoit,uit,tempi,yei,abari,        &
     &                 ptot,cs,etai,stot)                               
                  t9=tempi 
               endif 
            end if 
         enddo 
         write(*,*)'maximum burning steps exceeded! ' 
   92    continue 
!                                                                       
!--update particle quantities (rescale composition)                     
!                                                                       
         call epsb(yinput,ya,nel,enrg) 
         enrg=enrg/uergg 
         u(icell)=u(icell) + enrg 
         temp(icell)=t9 
         scale=1./(1.+err) 
         xpn=1.d0 
         do i=1,nel 
            ycc(icell,i+3)=real(scale*ya(i)) 
            xpn=xpn-scale*ya(i)*amas(i+3) 
         enddo 
         xpn=max(0.d0,xpn) 
         xp(icell)=ye(icell)*xpn 
         xn(icell)=(1.d0-ye(icell))*xpn 
         stepnuc=stepnuc + enrg*deltam(icell) 
   95    continue 
      enddo 
      totnuc=totnuc + stepnuc 
      write(*,*)'total nuc. en. generated so far : ',totnuc 
      write(*,*)'total nuc. en. generated this dt: ',stepnuc 
      write(80,*)time,stepnuc 
      call flush(80) 
!                                                                       
      return 
      END                                           
      double precision function accur(y,n) 
      implicit double precision (a-h,o-z) 
      dimension y(15),a(14) 
!---------------------------------------------------------------        
! program checks accuracy of mass conservation sum (a(i)*y(i)=1)        
!---------------------------------------------------------------        
      data a/4.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,            &
     &       56.,60./                                                   
      summ=0.0 
      do 10 i=1,n 
      summ=summ+a(i)*y(i) 
   10 continue 
      accur=summ - 1.0 
      return 
      END                                           
      subroutine derivn(y,rob,ff,rrat,rlam) 
!-----------------------------------------------------------------      
! calculates time derivatives of abundances i.e.                        
! ff(i) is system of differential equations                             
!-----------------------------------------------------------------      
      implicit double precision (a-h,o-z) 
      dimension rrat(17),rlam(15),y(15),ff(14) 
      do 10 i=1,14 
      ff(i)=0.0 
   10 continue 
      y(15)=0.0 
      do 20 i=2,14 
      ff(1)=ff(1)-rob*rrat(i)*y(1)*y(i)+rlam(i+1)*y(i+1) 
   20 continue 
      ff(1)=ff(1)-0.5*rrat(1)*rob**2*y(1)**3                            &
     &      +3.0*rlam(2)*y(2)                                           &
     &      +0.5*rrat(15)*rob*y(2)**2                                   &
     &      +rrat(16)*rob*y(2)*y(3)                                     &
     &      +0.5*rrat(17)*rob*y(3)**2                                   
      ff(2)=1./6.*rrat(1)*rob**2*y(1)**3+rlam(3)*y(3)                   &
     &      -rrat(2)*rob*y(1)*y(2)-rlam(2)*y(2)                         &
     &      -rrat(15)*rob*y(2)**2                                       &
     &      -rrat(16)*rob*y(2)*y(3)                                     
      do 30 i=3,14 
      ff(i)=rrat(i-1)*rob*y(1)*y(i-1)+rlam(i+1)*y(i+1)                  &
     &      -rrat(i)*rob*y(1)*y(i)-rlam(i)*y(i)                         
   30 continue 
      ff(3)=ff(3)-rrat(16)*rob*y(2)*y(3)-rrat(17)*rob*y(3)**2 
      ff(4)=ff(4)+0.5*rrat(15)*rob*y(2)**2 
      ff(5)=ff(5)+rrat(16)*rob*y(2)*y(3) 
      ff(6)=ff(6)+0.5*rrat(17)*rob*y(3)**2 
      return 
      END                                           
      subroutine epsb(yold,y,n,enrg) 
!--------------------------------------------------------------         
! change of total nuclear binding energy/mass                           
!       units erg/g                                                     
!................................................................       
      implicit double precision (a-h,o-z) 
      dimension yold(15),y(15),amex(14) 
      common /mass/ amex 
      data amex /2.425, 0.0, -4.737, -7.046, -13.933, -21.492,          &
     &           -26.016, -30.231, -34.845, -37.549, -42.818,           &
     &           -48.331, -53.902, -54.185/                             
      enrg=0.0 
      do 10 i=1,n 
   10 enrg=enrg-(y(i)-yold(i))*amex(i) 
      enrg=9.616221376d17*enrg 
      if(abs(enrg).lt.1.d4) enrg=0. 
      return 
      END                                           
      subroutine genpar(rho,t9,a,z,y,ye,ar,ar11,coefg,coeft) 
!*************************************************************          
!                                                            *          
!  this subroutine computes the screening parameters for the *          
!  Thielemann network.                                       *          
!                                                            *          
!*************************************************************          
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      dimension y(15),a(14),z(14),ar(14) 
      ye=0.0 
      do 20 i=1,14 
      ye=ye+z(i)*y(i) 
   20 continue 
      coefa=7.345889d-9*(rho*ye)**(0.33333333333) 
      do 30 i=1,14 
      ar(i)=coefa*z(i)**(0.3333333333) 
   30 continue 
      ar11=coefa*4.**(0.3333333333) 
      coefg=1.67100d-12/t9 
      coeft=3.58039d4*coefg**(0.3333333333333) 
      return 
      END                                           
      subroutine lequb (a,b,n) 
!-------------------------------------------------------------          
!     program solves system of linear equations a*x=b                   
!     solution x in b                                                   
!-------------------------------------------------------------          
      implicit double precision (a-h,o-z) 
      dimension a(  n,  n),B(  n) 
      n1=n-1 
!     search for maximum element in each row and divide                 
      do 19 i=1,n 
      r=dabs(a(i,1)) 
      do 16 j=2,n 
      c=dabs(a(i,j)) 
   16 if(c.gt.r) r=c 
      do 17 j=1,n 
   17 a(i,j)=a(i,j)/r 
  824 b(i)=b(i)/r 
   19 continue 
      do 9 j=1,n1 
      l=j+1 
!     now pivoting                                                      
    6 do 9 i=l,n 
      r=-a(i,j)/a(j,j) 
      if (r.eq.0.0) goto 9 
      do 7 k=l,n 
    7 a(i,k)=a(i,k)+r*a(j,k) 
  825 b(i)=b(i)+r*b(j) 
    9 continue 
!     matrix now in triangular from => solution                         
  826 b(n)=b(n)/a(n,n) 
      do 13 l=1,n1 
      i=n-l 
      r=0.d0 
      imax=i+1 
      do 12 j=imax,n 
      jj=i+n+1-j 
   12 r=r+a(i,jj)*b(jj) 
   13 b(i)=(b(i)-r)/a(i,i) 
   15 return 
      END                                           
      subroutine matr(y,rob,rrat,rlam,dt,g) 
!--------------------------------------------------------------         
! g(i,j) matrix for multi-dimensional newton-raphson                    
! g(i,j)=dff(i)/dy(j)-delta(i,j)/deltat                                 
!--------------------------------------------------------------         
      implicit double precision (a-h,o-z) 
      dimension rrat(17),rlam(15),y(15),g(14,14) 
      do 10 i=1,14 
      do 10 j=1,14 
      g(i,j)=0.0 
   10 continue 
      do 50 i=2,14 
      term=rrat(i)*rob 
      g(1,i)=g(1,i)-term*y(1)+rlam(i) 
      g(1,1)=g(1,1)-term*y(i) 
   50 continue 
      g(1,2)=g(1,2)+2.0*rlam(2)+rrat(15)*rob*y(2)+rrat(16)*rob*y(3) 
      g(1,3)=g(1,3)+rrat(16)*rob*y(2)+rrat(17)*rob*y(3) 
      g(1,1)=g(1,1)-1.5*rrat(1)*rob**2*y(1)**2 
      g(2,1)=0.5*rrat(1)*rob**2*y(1)**2-rrat(2)*rob*y(2) 
      g(2,2)=-rrat(2)*rob*y(1)-rlam(2)-2.0*rrat(15)*rob*y(2)            &
     &       -rrat(16)*rob*y(3)                                         
      g(2,3)=rlam(3)-rrat(16)*rob*y(2) 
      do 60 i=3,14 
      g(i,i-1)=g(i,i-1)+rrat(i-1)*rob*y(1) 
      g(i,i)=g(i,i)-rrat(i)*rob*y(1)-rlam(i) 
      g(i,1)=rrat(i-1)*rob*y(i-1)-rrat(i)*rob*y(i) 
   60 continue 
      do 70 i=3,13 
      g(i,i+1)=g(i,i+1)+rlam(i+1) 
   70 continue 
      g(3,2)=g(3,2)-rrat(16)*rob*y(3) 
      g(3,3)=g(3,3)-rrat(16)*rob*y(2)-2.0*rrat(17)*rob*y(3) 
      g(4,2)=g(4,2)+rrat(15)*rob*y(2) 
      g(5,2)=g(5,2)+rrat(16)*rob*y(3) 
      g(5,3)=g(5,3)+rrat(16)*rob*y(2) 
      g(6,3)=g(6,3)+rrat(17)*rob*y(3) 
      do 80 i=1,14 
      g(i,i)=g(i,i)-1.0/dt 
   80 continue 
      return 
      END                                           
      subroutine newab(yold,y,rob,rrat,rlam,ff,dt,err,eps) 
      implicit double precision (a-h,o-z) 
      dimension yold(15),y(15),rrat(17),rlam(15),ff(14) 
      dimension rs(14),g(14,14) 
      data maxit/100/ 
!-------------------------------------------------------------          
! this subroutine computes the new abundances y(i) at t+dt              
! yold and dt have to be given                                          
! rob (density in g/cm**3), rrat(i), rlam(i), ff(i) are                 
! input ; the rates and derivatives have to be provided                 
! they are calculated with rates and deriv for conditions               
! (t9,rob,y(i)) at the beginning of the time step                       
! err is the accuracy of mass conservation, eps is the allowed          
! error in the newton-raphson iteration                                 
! ---- you can play with the size of err and eps ----------             
!-------------------------------------------------------------          
      dt1=1./dt 
      do 30 i=1,14 
      y(i)=yold(i) 
   30 continue 
      knew=0 
   35 call matr(y,rob,rrat,rlam,dt,g) 
      knew=knew + 1 
      do 90 i=1,14 
      rs(i)=(y(i)-yold(i))*dt1 - ff(i) 
   90 continue 
      ers=0.0 
      call lequb(g,rs,14) 
      do 100 i=1,14 
!      write (51,*) i,ff(i),y(i),yold(i)                                
      y(i)=y(i)+rs(i) 
!                                                                       
!--if abundance negative set equal to zero                              
!  (this has to be commented out if four lines down are                 
!                                                                       
      if(y(i).le.0.0d0)then 
         y(i)=0.0d0 
         go to 100 
      else 
         ersi=dabs(rs(i)/y(i)) 
         ers=dmax1(ers,ersi) 
      end if 
  100 continue 
      err=accur(y,14) 
!      testsum=0.                                                       
!      do j=1,14                                                        
!         testsum=testsum+y(j)                                          
!         print *, y(j),testsum                                         
!      end do                                                           
!                                                                       
!--iterate on solution until satisfy error tolerance                    
!                                                                       
      if(dabs(err).lt.eps) goto 110 
      if(ers.lt.eps) goto 110 
      call derivn(y,rob,ff,rrat,rlam) 
      if(knew.gt.maxit)then 
         write(*,*)' error: knew gt maxit' 
         stop 
      end if 
      goto 35 
  110 return 
      END                                           
      subroutine rates(icell,t9,pf,pb,rrat,rlam,ye,coefg,coeft,a,       &
     &                 ar,ar11,z)                                       
      implicit double precision (a-h,o-z) 
      dimension pf(7,17),pb(7,14),rrat(17),rlam(15),a(14),              &
     &          ar(14),z(14), ptemp(7)                                  
!--------------------------------------------------------               
! calculates rates                                                      
! rrat(i) forward rates (captures 1-14) 15:CC 16:CO 17:OO               
! pf(7,i) coefficients                                                  
! rlam(i) backward rates (1-14 photo-disintegrations)                   
! pb(7,i) coefficients                                                  
!--------------------------------------------------------               
      t91=t9**(-1) 
      t92=t9**(-1./3.) 
      t93=t9**(1./3.) 
      t94=t9 
      t95=t9**(5./3.) 
      t96=dlog(t9) 
      rrat(14)=0.0 
      rlam(15)=0.0 
      rlam(1)=0.0 
      do 10 i=2,13 
         do it=1,7 
            ptemp(it)=pb(it,i) 
         enddo 
         rlam(i)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
         do it=1,7 
            ptemp(it)=pf(it,i) 
         enddo 
         rrat(i)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
         rrat(i)=rrat(i)*scrng(z(1),z(i),a(1),a(i),ar(1),ar(i),ye,      &
     &                       coefg,coeft)                               
   10 continue 
      do it=1,7 
         ptemp(it)=pf(it,1) 
      enddo 
      rrat(1)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
      rrat(1)=rrat(1)*scrng(z(1),z(1),a(1),a(1),ar(1),ar(1),ye,         &
     &                       coefg,coeft)                               
      rrat(1)=rrat(1)*scrng(z(1),4.d0,a(1),8.d0,ar(1),ar11,             &
     &                       ye,coefg,coeft)                            
      do it=1,7 
         ptemp(it)=pf(it,15) 
      enddo 
      rrat(15)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
      rrat(15)=rrat(15)*scrng(z(2),z(2),a(2),a(2),ar(2),ar(2),          &
     &                            ye,coefg,coeft)                       
      do it=1,7 
         ptemp(it)=pf(it,16) 
      enddo 
      rrat(16)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
      rrat(16)=rrat(16)*scrng(z(2),z(3),a(2),a(3),ar(2),ar(3),          &
     &                            ye,coefg,coeft)                       
      do it=1,7 
         ptemp(it)=pf(it,17) 
      enddo 
      rrat(17)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
      rrat(17)=rrat(17)*scrng(z(3),z(3),a(3),a(3),ar(3),ar(3),          &
     &                            ye,coefg,coeft)                       
      do it=1,7 
         ptemp(it)=pb(it,14) 
      enddo 
      rlam(14)=rk(ptemp,t91,t92,t93,t94,t95,t96) 
      return 
      END                                           
      subroutine rrate(pf,pb,iprint,idisk2) 
      implicit double precision (a-h,o-z) 
      dimension pf(7,17),pb(7,14) 
!------------------------------------------------------------           
! program reads coefficients of rates                                   
! pf forward rates, pb backward rates (=photodisintegrations)           
!------------------------------------------------------------           
!                                                                       
      character*72 string 
      do 10 i=1,17 
      if(i.eq.14) goto 10 
      read(idisk2,1000)string 
!      write(*,1000)string                                              
      read(idisk2,2000) (pf(j,i),j=1,7) 
   10 continue 
      do 20 i=2,14 
      read(idisk2,1000)string 
!      write(*,1000)string                                              
      read(idisk2,2000) (pb(j,i),j=1,7) 
   20 continue 
 1000 format(a72) 
 2000 format(4d13.6) 
      return 
      END                                           
      double precision function scrng(z1,z2,a1,a2,ar1,                  &
     &                 ar2,ye,coefg,coeft)                              
      implicit double precision (a-h,o-z) 
      gamma=coefg*z1*z2*2./(ar1+ar2) 
      tau=coeft*(a1*a2/(a1+a2)*z1**2*z2**2)**(0.33333333333) 
      scrng=dexp(1.25*gamma-0.0975*tau*(3.*gamma/tau)**2) 
      return 
      END                                           
      function rk(p,t91,t92,t93,t94,t95,t96) 
      implicit double precision (a-h,o-z) 
      dimension p(7) 
      rk=p(1)+p(2)*t91+p(3)*t92+p(4)*t93+p(5)*t94+p(6)*t95              &
     &   +p(7)*t96                                                      
      rk=dexp(rk) 
      return 
      END