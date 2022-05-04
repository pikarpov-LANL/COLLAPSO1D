      program lahyc 
!*******************************************************                
!                                                      *                
!  This is a Lagrangian 1D hydrodynamics code.         *                
!  adds particles                                      *                
!  this is the main driver of the integration.         *                
!                                                      *                
!*******************************************************                
!
      !use torch_ftn
      !use iso_fortran_env
      implicit double precision (a-h,o-z) 
!                                                                       
!--ntstep counts the number of timesteps                                
      integer jtrape,jtrapb,jtrapx,ntstep 
!                                                                       
      parameter (idim=4000) 
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
      logical trapnue, trapnueb, trapnux, print_nuloss 
      common /trap / trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /typef/ iextf, ieos 
      common /units/ umass, udist, udens, utime, uergg, uergcc 
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg 
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1 
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu 
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
!--nstep is a counter for dump fiels                                    
!--ntstep is the number of timesteps                                    
      nstep = 0 
      ntstep = 1 
                                                                        
      do while (time.lt.tmax) 
!42   continue                                                          
!                                                                       
        nstep = nstep + 1 
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
        print*,' ********' 
        print*,' savetime = ',time 
        print*,' ********' 
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
      call printout(lu) 
      stop 
      END                                           
      subroutine hydro(time,ncell,x,v,                                  &
     &           u,rho,ye,f,du,dye,q,fmix,                              &
     &           ynue,ynueb,ynux,dynue,dynueb,dynux,                    &
     &           unue,unueb,unux,dunue,dunueb,dunux,                    &
     &           print_nuloss)                                          
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
      integer jtrape,jtrapb,jtrapx 
      parameter (idim=4000) 
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
      !print *, '[time] ',time,ieos                                     
      if(ieos.eq.1)call eospg(ncell,rho,u) 
      if(ieos.eq.2)call eospgr(ncell,rho,u) 
      if(ieos.eq.3.or.ieos.eq.4)call eos3(ncell,rho,u,ye) 
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
         if (ynue(k).lt.tiny) print *,'yne<0, k=',k,ynue(k),            &
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
!--turbulence contribution to pressure via ML 
      !call testml
      call turbpress(ncell,rho) 
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
      call turbulence(ncell,x,f,q,v,rho,fmix) 
!                                                                       
!--compute q values                                                     
!                                                                       
      call artvis(ncell,x,rho,v,q) 
!                                                                       
!--compute forces on the particles                                      
!                                                                       
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
!                                                                       
      dimension x(0:idim),rho(idim) 
      common /carac/ deltam(idim), abar(idim) 
!                                                                       
      data pi4/12.566371/ 
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
      parameter (idim=4000) 
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
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
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
         vturbk=dsqrt(vturb2(kp05)) 
         akp1=pi4*x(k1)*x(k1) 
         akp=pi4*x(k)*x(k) 
         pdv=(pr(kp05)+rho(kp05)*vturb2(kp05))*                         &
     &        (akp1*v(k1)-akp*v(k))/deltam(kp05)                        
!                                                                       
         dupp=0.0 
         if (temp(kp05).lt.6.and.rho(kp05).lt.1000.) then 
            dupp=9.96d5*temp(kp05)**9/rho(kp05) 
         end if 
         if (ifleos(kp05).ne.3) then 
!-- we subtract the energy added to the neutrino field                  
! no shock heating                                                      
            du(kp05)=-pdv+0.5*dq(kp05)-dunu(kp05)-dupp+                 &
     &           rho(kp05)*vturbk**3/dmix(kp05)                         
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension x(0:idim),f(0:idim),v(0:idim) 
      dimension q(idim),rho(idim) 
!                                                                       
      common /eosnu/ prnu(idim1) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /carac/ deltam(idim), abar(idim) 
      common /damping/ damp, dcell 
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
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
         if (rho(km05+1).gt.0.1.and.km05.lt.1.and.                      &
     &        1.5d0*rho(km05)*xk3.lt.rho(kp05)*xkp3) then               
            pressp=(pr(kp05)+prnu(kp05)) 
            print *, 'increased force' 
         else 
            pressp = pr(kp05) + prnu(kp05) 
         end if 
         pressm = pr(km05) + prnu(km05) 
         gradp=ak*(pressp - pressm) 
         gradpt=ak*(rho(kp05)*vturb2(kp05)-                             &
     &        rho(km05)*vturb2(km05))                                   
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
      parameter (idim=4000) 
!                                                                       
      dimension x(0:idim),f(0:idim) 
      dimension gpot(idim),deltam(idim) 
      dimension xmi(0:idim) 
!                                                                       
      common /fturb/ geff(idim), fmix1(idim), fmix2(idim) 
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
         geff(k)=-gg*xmi(k)/r2 
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
      parameter (idim=4000) 
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
      parameter(idim=4000) 
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
      parameter(idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      parameter(idim=4000) 
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
      parameter (idim=4000) 
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
      parameter(idim=4000) 
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
      parameter (idim=4000) 
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
      f=ufoe/utime 
      if (print_nuloss.eqv..true.) then 
        write(*,510)'[nue loss, foes/s (MeV)]',rlumnue*f,               &
     &  '               ',enue                                          
  510   format(A,1p,E10.3,A,E10.3) 
        !print*,'nueb loss:',rlumnueb*f,' foes/s at <E> (MeV)',enueb    
        !print*,' nux loss:',rlumnux*f,' foes/s at <E> (MeV)',enux      
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
      parameter(idim=4000) 
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
      parameter(idim=4000) 
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
      subroutine turbpress(ncell,rho) 
!*******************************************************                
!                                                                       
! This subroutine passes the grid to a trained PyTorch model            
! that predicts pressure due to turbulence                              
!                                                                       
!******************************************************                 
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      parameter(idim=4000) 
      parameter (idim1=idim+1) 
!                                                                       
      dimension rho(idim) 
      dimension unue(idim),unueb(idim),unux(idim) 
!                                                                       
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim) 
      common /eosnu/ prnu(idim1) 
      common /turb/ prturb(idim1) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
!                                                                       
      pturb = 1 
      !do i=1,ncell                                                     
      !  prturb(i)=                                                     
      !enddo                                                            
      !print *,'nupress: maximum prnu/pr,etanuemx',ratmax,etanuemx      
!                                                                       
      return 
      END                                           
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
      parameter(idim=4000) 
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
      parameter(idim=4000) 
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
      parameter (idim=4000) 
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
      parameter (idim=4000) 
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
      rhoswe=1.e11/udens 
!      rhoswe=1.e11/udens                                               
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
      parameter(idim=4000) 
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
      parameter(idim=4000) 
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
!      print *,'rootemp: iteration did not converge'                    
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
!******************************************************                 
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
      if (ye.gt..5d0) print *, kcell,ye 
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
      subroutine readini 
!*************************************************************          
!                                                                       
!     reads initial conditions                                          
!                                                                       
!************************************************************           
!                                                                       
      implicit double precision (a-h,o-z) 
!                                                                       
      integer jtrape,jtrapb,jtrapx 
!                                                                       
      parameter (idim=4000) 
      parameter (idim1=idim+1) 
      parameter (iqn=17) 
      real ycc,yccave 
!                                                                       
      common /cc   / ycc(idim,iqn), yccave(iqn) 
      common /ener1/ dq(idim), dunu(idim) 
      common /numb/ ncell, ncell1 
      common /celle/ x(0:idim),v(0:idim),f(0:idim) 
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim) 
      common /nustuff/ ynue(idim),ynueb(idim),ynux(idim),               &
     &               unue(idim),unueb(idim),unux(idim)                  
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
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
!                                                                       
      character*9 filin,filout 
      data pi4/12.56637d0/ 
      gg=13.34 
      tacr=1.d2 
!                                                                       
!--read options                                                         
!                                                                       
      open(11,file='inlahyc') 
      read(11,10) filin 
      read(11,10) filout 
   10 format(A) 
      read(11,*) idump 
      read(11,*) dtime,tmax 
      read(11,*) cq,cl 
      read(11,*) iextf,ieos,dcore 
      read(11,*) ncell,delp,nups,damp,dcell 
      read(11,*) iflxlm, icvb, ufact, yefact 
      print *,'cq,cl',cq,cl 
      print *,'iextf,ieos',iextf,ieos 
!                                                                       
!--open binary file containing initial conditions                       
!                                                                       
                                                                        
      open(60,file=filin,form='unformatted') 
      open(61,file=filout,form='unformatted') 
!                                                                       
!--position pointer in binary file                                      
!                                                                       
      do i=1,idump-1 
         read(60) idummy 
      enddo 
!                                                                       
!--read data                                                            
!                                                                       
      nqn=17 
!                                                                       
      read(60) nc,t,xmcore,rb,ftrape,ftrapb,ftrapx,                     &
     &   (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),      &
     &      (u(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),          &
     &      (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),            &
     &      (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),           &
     &      (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),        &
     &      (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),        &
     &      (ufreez(i),i=1,nc),(pr(i),i=1,nc),(u2(i),i=1,nc),           &
     &      (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),              &
     &     (vturb2(i),i=1,nc),                                          &
     &     ((ycc(i,j),j=1,nqn),i=1,nc)                                  
!                                                                       
      time = t 
!                                                                       
      print *, nc, t, xmcore, rb, ftrape,ftrapb,ftrapx 
      do k=1,nc 
         write(44,*)k,temp(k),u(k),u2(k) 
      end do 
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
      do k=1,3 
         print *,u2(k),pr(k),ufreez(k) 
         print *, unue(k),unueb(k),unux(k) 
         print *, ynue(k),ynueb(k),ynux(k) 
         print *, xp(k), xn(k), ifleos(k), xmcore 
         do j=1,nqn 
            print *,ycc(k,j) 
         end do 
      end do 
      p1out=pr(ncell) 
      p2out=1.39*gg*deltam(ncell)/x(ncell)**4/pi4 
      rout=x(ncell) 
      print *, 'ncell = ', ncell 
      print *, 'xmcore = ', xmcore 
      print *, 'rho(1) = ', rho(1) 
      print *, x(ncell),pr(ncell) 
      print *, 'deltam(1) = ', deltam(1) 
      print *, 'time = ',time 
!      gamma=1.666666666666667                                          
      ncell1=ncell+1 
      if (ieos.eq.4) then 
         do i=1,ncell 
            write(71,*)i,u2(i),u(i),abar(i) 
            if (ifleos(i).eq.3) then 
               s=u2(i) 
               u2(i)=u(i) 
               u(i)=s 
               write (71,*) i,u2(i),u(i),abar(i) 
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
!                                                                       
      parameter (idim=4000) 
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
      common /state/ xp(idim), xn(idim), eta(idim), ifleos(idim) 
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax 
      common /freez/ ufreez(idim) 
      common /core / dcore, xmcore 
      common /typef/ iextf, ieos 
      common /bstuf/ rb, dumrb, f1rb, f2rb 
      common /tempe/ temp(idim) 
      common /cgas / gamma 
      common /timej / time, dt 
      logical trapnue, trapnueb, trapnux 
      common /trap/ trapnue(idim), trapnueb(idim), trapnux(idim) 
      common /ftrap/ ftrape,ftrapb,ftrapx 
      common /jtrap/ jtrape,jtrapb,jtrapx 
      common /carac/ deltam(idim), abar(idim) 
      common /ener2/ tkin, tterm 
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
      logical te(idim), teb(idim), tx(idim) 
      dimension uint(idim), s(idim) 
      equivalence(trapnue,te) 
      equivalence(trapnueb,teb) 
      equivalence(trapnux,tx) 
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
!                                                                       
!      do i=1,1                                                         
!         write (12,*) i, nc,t,gamma,tkin,tterm                         
!      end do                                                           
!                                                                       
      write(lu)nc,t,xmcore,rb,ftrape,ftrapb,ftrapx,                     &
     &   (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),      &
     &      (uint(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),       &
     &      (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),            &
     &      (xp(i),i=1,nc),(xn(i),i=1,nc),(ifleos(i),i=1,nc),           &
     &      (ynue(i),i=1,nc),(ynueb(i),i=1,nc),(ynux(i),i=1,nc),        &
     &      (unue(i),i=1,nc),(unueb(i),i=1,nc),(unux(i),i=1,nc),        &
     &      (ufreez(i),i=1,nc),(pr(i),i=1,nc),(s(i),i=1,nc),            &
     &     (te(i),i=1,nc),(teb(i),i=1,nc),(tx(i),i=1,nc),               &
     &     (vturb2(i),i=1,nc),                                          &
     &     ((ycc(i,j),j=1,nqn),i=1,nc)                                  
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
      integer jtrape,jtrapb,jtrapx, ntstep 
!                                                                       
      parameter (idim=4000) 
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
      common /fturb/ geff(idim), fmix1(idim), fmix2(idim) 
      common /turb/ vturb2(idim),dmix(idim),alpha(4),bvf(idim) 
      common /dturb/ dumvt2(idim) 
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
!                                                                       
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
      if(ifirst) then 
         ifirst=.false. 
         dumx(0)=x(0) 
         do i=1,ncell 
            steps(i)=dt0/1.d8 
            istep(i)=nint(dtime/steps(i)) 
            t0(i)=time 
            tempold(i)=temp(i) 
            rhold(i)=rho(i) 
            yeold(i)=ye(i) 
            xnold(i)=xn(i) 
            xpold(i)=xp(i) 
         enddo 
         print *, 'calling hydro first' 
         print_nuloss=.false. 
         call hydro(time,ncell,x,v,                                     &
     &           u,rho,ye,f1v,f1u,f1ye,q,fmix1,                         &
     &           ynue,ynueb,ynux,f1ynue,f1ynueb,f1ynux,                 &
     &           unue,unueb,unux,f1unue,f1unueb,f1unux,                 &
     &           print_nuloss)                                          
      end if 
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
            dumvt2(i)=vturb2(i) + dtf11*fmix1(i) 
            reset(i)=.false. 
            if (dumye(i).lt.0.02) then 
               dumye(i)=.02 
            end if 
         enddo 
!                                                                       
!--get forces at half the time step, using predictions                  
!                                                                       
         thalf=time + f11*steps(1) 
         print_nuloss=.true. 
         call hydro(thalf,ncell,dumx,dumv,                              &
     &         dumu,rho,dumye,f2v,f2u,f2ye,q,fmix2,                     &
     &         dumynue,dumynueb,dumynux,f2ynue,f2ynueb,f2ynux,          &
     &         dumunue,dumunueb,dumunux,f2unue,f2unueb,f2unux,          &
     &         print_nuloss)                                            
!                                                                       
!--advance all the gas particles                                        
!                                                                       
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
            dumvt2(i)=vturb2(i) + dtf21*fmix1(i) + dtf22*fmix2(i) 
            if (dumye(i).lt.0.02) then 
               dumye(i)=.02 
            end if 
         enddo 
!                                                                       
!     set saturation const=0                                            
         satc=0 
!                                                                       
!--get forces at end of time step                                       
!                                                                       
         tfull=time + steps(1) 
         print_nuloss=.false. 
         call hydro(tfull,ncell,dumx,dumv,                              &
     &         dumu,rho,dumye,f2v,f2u,f2ye,q,fmix2,                     &
     &         dumynue,dumynueb,dumynux,f2ynue,f2ynueb,f2ynux,          &
     &         dumunue,dumunueb,dumunux,f2unue,f2unueb,f2unux,          &
     &         print_nuloss)                                            
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
               else 
                  if (ieos.eq.3) then 
                     usi=u(i)+8.8*umevnuc 
                  else 
!--variable of state is entropy                                         
                     usi=u(i)*temp(i) 
                  endif 
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
            vturb2(i)=dumvt2(i) 
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
     &              u,rho,ye,f1v,f1u,f1ye,q,fmix1,                      &
     &              ynue,ynueb,ynux,f1ynue,f1ynueb,f1ynux,              &
     &              unue,unueb,unux,f1unue,f1unueb,f1unux,              &
     &              print_nuloss)                                       
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
         call eosflg(ncell,rho,ye,u,f1ye,f1u) 
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
            write(*,500)'[    time/tmax, dt     ]',                     &
     &           time,'/',tmax,steps(1)                                 
  500       format(A,1p,E10.3,A,E9.3,E15.3) 
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
      parameter (idim=4000) 
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



!---------ocean--------

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
! *** Initialize variables to avoid fp exceptions                       
      g21=0.0 
      g2 = 0.0 
      g41=0.0 
      g4=0.0 
!                                                                       
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
!     Avoid floating point underflow                                    
      if(pc*z4.lt.-100) then 
         z1 = 0.0 
      else 
         z1=exp(pc*z4) 
      endif 
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

!------nse4c------



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
     &                   xa,xh,yeh,zbar,abar)                           
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
     &                   xa,xh,yeh,zbar,abar)                           
!         write(*,999) t9fnl,t9,rho,ye,yp,yn,kit                        
          If(kit.ge.kmtmp) Then 
             write(*,999) t9fnl,t9,rho,ye,yp,yn,kit 
             print*,'broken code' 
!            stop                                                       
          Endif 
      Endif 
      Return 
  999 Format(6(1pd11.4),1x,I3) 
                                                                        
      END                                           
                                                                        
                                                                        
      subroutine nsesolv (ider,t9,rho,ye,yp,yn,kit,kmax,ubind,          &
     &                   xa,xh,yeh,zbar,abar)                           
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
              yt=dexp(max(dlgalp(n)+zpt*dlgyp+znt*dlgyn+h(iz), -50.)) 
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
          ytst2=0.0d0                                                   
          benuc=0.d0                                                    
          dt9=1.d0/t9                                                   
!--do loop in two parts for speed                                       
          DO m=1,3                                                      
            iz=intz(m)                                                  
            zpm=zp(m)                                                   
            zam=za(m)                                                   
            bim=bi(m)                                                   
            ym=dexp(max(dlgalp(m)+dlgyp*zpm+dlgyn*zn(m)+h(iz), -50.))   
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
          ENDDO                                                         
!--alpha particle mass fraction (last xm computed above)                
          xa=xm                                                         
          xh=0.0d0                                                      
          yeh=0.0d0                                                     
          DO m=4,nmax                                                   
            iz=intz(m)                                                  
            zpm=zp(m)                                                   
            zam=za(m)                                                   
            bim=bi(m)                                                   
            ym=dexp(max(dlgalp(m)+dlgyp*zpm+dlgyn*zn(m)+h(iz), -50.))   
!           y(m)=ym                                                     
            xm=zam*ym                                                   
!--heavies mass fraction                                                
            xh=xh+xm                                                    
!           yeh=yeh+xm*zpm/zam                                          
            yeh=yeh+zpm*ym                                              
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
          ENDDO                                                         
!--normalize yeh                                                        
          if (xh.lt.1d-10) then                                         
            xh=0.0d0                                                    
            yeh=0.0d0                                                   
          else                                                          
            yeh=yeh/xh                                                  
          endif                                                         
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
      dHtdT9=0.0d0 
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

!--------sleos--------

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
          H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4 
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
        stop 
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
   11    U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))            &
     & /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P                    
      Y2(N)=0. 
      DO 12 K=N-1,1,-1 
   12    Y2(K)=Y2(K)*Y2(K+1)+U(K) 
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
    1 IF(KHI-KLO.EQ.1) GOTO 2 
      K=(KHI+KLO)/2 
      IF(XA(K).GT.X)THEN 
      KHI=K 
      ELSE 
      KLO=K 
      ENDIF 
      GOTO 1 
    2 H=XA(KHI)-XA(KLO) 
      IF(H.EQ.0.) PAUSE 'BAD XA INPUT. ' 
!c  Compute spline fit.                                                 
      A=(XA(KHI)-X)/H 
      B=(X-XA(KLO))/H 
      Y=A*YA(KLO)+B*YA(KHI)+                                            &
     & ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H**2/6.                    
!      write(*,5)klo,khi,x,xa(klo),xa(khi),ya(klo),ya(khi)              
!     > ,y2a(klo),y2a(khi),y                                            
    5   format(2i3,1p8e9.2) 
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
                          