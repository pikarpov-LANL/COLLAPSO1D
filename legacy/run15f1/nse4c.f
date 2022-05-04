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
c     write(*,*) 'testt',testt
      If(dabs(testt).gt.nsetol) Then
          call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmax,ubind,
     &                   xa,xh,yeh,zbar,abar)
c         write(*,999) t9fnl,t9,rho,ye,yp,yn,kit
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
          call nsesolv(ider,t9,rho,ye,yp,yn,kit,kmtmp,ubind,
     &                   xa,xh,yeh,zbar,abar)
c         write(*,999) t9fnl,t9,rho,ye,yp,yn,kit
          If(kit.ge.kmtmp) Then
             write(*,999) t9fnl,t9,rho,ye,yp,yn,kit
             print*,'broken code'
c            stop
          Endif
      Endif
      Return
  999 Format(6(1pd11.4),1x,I3)

      End


      subroutine nsesolv (ider,t9,rho,ye,yp,yn,kit,kmax,ubind,
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
c     write(*,*) 'NSEsolv',t9,rho,ye
c--Raph's routine to find abar, zbar, and nuclear binding energy
c  ---------------------------------------------------------------------------
c  This portion is designed to calculate the Nuclear Statistical equilibrium
c  distribution of material for a given temperature and density.  Analytic
c  expressions exist which give solutions for Y(A,Z) as a product of separable
c  functions of Yp, Yn, and alpha, where alpha contains all of the temperature,
c  density, and nuclear parameter dependence.  Yp and Yn are the free proton
c  and neutron abundances.  Implicit solution of equations for Ye and mass
c  conservation, using the analytic expressions to convert Y(A,Z) to Yp and Yn,
c  yields values for Yp and Yn and hence the entire distribution.
c  We will employ a two dimensional Newton method, which, though it converges
c  extremely rapidly if it is near the solution, is prone to severe 
c  convergence problems.  For this reason the subroutine uses the values of
c  Yp and Yn from the previous step.  The subroutine also reports if it 
c  converges.  In the event that it fails to converge, the correct solution 
c  requires a more detailed iteration in T, rho and Ye.
c  ---------------------------------------------------------------------------
      dlgyp=dlog(yp)
      dlgyn=dlog(yn)
      Do 100 i=1,nmax
        If(zp(i).gt.1.d0) Then
            intz(i)=int(zp(i))
        Else
            intz(i)=1
        Endif
  100 Continue
c--Calculate the thermodynamic dependence coefficients
      call alphcalc(nmax,t9,rho)
c     call nsenorm(y,za,nmax,nmax,zp,zn,ye)
c--Calculate the screening factors
      call nsescreen(T9,rho,ye,h,dhdp,dhdn,dhdT9,iscrn)
      k=1
      testk=1.d0+tol
c  ---------------------------------------------------------------------------
c  For each iterative T and rho it is necessary to first calculate alph(n), 
c  which contains the dependence on temperature, density and nuclear parameters
c  In the case of screening (iscrn=2) the alph(n) are then modified by the 
c  appropriate screening factor, the product of the screening factors for 
c  the series of (p,gamma) reactions necessary to build up the required Z.
c  Once alph(n) has been calculated, the next step is to calculate f, the 
c  implicit formula for Ye, g, the implicit formula for mass conservation and
c  dfp, dfn, dgn, and dgp, the derivatives of f and g with respect to yp and yn.
c  These 6 values, calculated using logaritmic interim steps to avoid overflow,
c  and the analytic solution of the matrix eqn, yield delyp and delyn. The 
c  new values of yp,yn are yp+delyp, yn+delyn. The implicit solution continues
c  until either kmax iterations are run, or the change in solution is less
c  than some minimum tolerance, tol.
c  ---------------------------------------------------------------------------
  20  If(k.lt.kmax) THEN
          IF(testk.gt.tol) THEN
            f=-ye
            g=-1.d0
            dfp=0.d0
            dfn=0.d0
            dgp=0.d0
            dgn=0.d0
c           dfp2=0.d0
c           dfn2=0.d0
c           dgp2=0.d0
c           dgn2=0.d0
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
c                 write(*,*) 'y explosion'
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
c             dfp2=dfp2+zpt*dhdp(iz)*yt
              dfn=dfn+zpt*zndyn
c             dfn2=dfn2+zpt*dhdn(iz)*yt
              dgp=dgp+zat*zpdyp
c             dgp2=dgp2+zat*dhdp(iz)*yt
              dgn=dgn+zat*zndyn
c             dgn2=dgn2+zat*dhdn(iz)*yt
  30        CONTINUE
c           write(*,*) 'K',k,f,g
c           dfn=dfn+dfn2
c           dfp=dfp+dfp2
c           dgn=dgn+dgn2
c           dgp=dgp+dgp2
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
c           write(10,4000)k,yp,yn
c           write(10,4300)f,g
c           write(15,4200)dfp,dfn,dgp,dgn
c           write(15,4250)dfp2,dfn2,dgp2,dgn2
c           write(10,4100)delyp,delyn,det
            k=k+1
            testzy=f
            testay=g
            testyp=delyp/yp
            testyn=delyn/yn
            testk=sqrt(testyp*testyp+testyn*testyn+f*f+g*g)
            GOTO 20
          ENDIF
c--NSE converges
          kit=k
c         write(10,4000)k,yp,yn
c         write(10,3100)T9r,rhor,ye
      ELSE
c--NSE fails to converge in Kmax steps
c         write(10,3200)kmax,T9r,rhor 
          yp=yp0
          yn=yn0
          kit=k
      ENDIF
c  -------------------------------------------------------------------------
c  Having completed the loop for T,rho (in kit < kmax iterations) or 
c  discovered that the solution will not converge (making kit >= kmax 
c  depending on how it fails) the subroutine returns c  for the next set 
c  of T, rho and Ye.  If ider >1 the subroutine calculates several moments 
c  of the distribution, including the average A and Z and the binding energy. 
c  -------------------------------------------------------------------------
      If(ider.ge.1) Then
          dkT2=1.d0/(8.6174d-2*t9*t9)
          atst=0.d0
          ztst=0.d0
          ytst=0.d0
          ytst2=0.0d0
          benuc=0.d0
          dt9=1.d0/t9
c--do loop in two parts for speed
          DO m=1,3
            iz=intz(m)
            zpm=zp(m)
            zam=za(m)
            bim=bi(m)
            ym=dexp(max(dlgalp(m)+dlgyp*zpm+dlgyn*zn(m)+h(iz), -50.))
c           y(m)=ym
            xm=zam*ym
c           x(m)=xm
            zpym=zpm*ym
c           at(iz)=at(iz)+xm
c           zt(iz)=zt(iz)+zpym
c
c--abar, zbar, benuc, are the average baryon number, proton number 
c  number, and binding energy per nucleon (in Mev/nuc).
c
            ytst=ytst+ym
            atst=atst+xm
            ztst=ztst+zpym
            benuc=benuc+bim*ym
          ENDDO
c--alpha particle mass fraction (last xm computed above)
          xa=xm
          xh=0.0d0
          yeh=0.0d0
          DO m=4,nmax
            iz=intz(m)
            zpm=zp(m)
            zam=za(m)
            bim=bi(m)
            ym=dexp(max(dlgalp(m)+dlgyp*zpm+dlgyn*zn(m)+h(iz), -50.))
c           y(m)=ym
            xm=zam*ym
c--heavies mass fraction
            xh=xh+xm
c           yeh=yeh+xm*zpm/zam
            yeh=yeh+zpm*ym
c           x(m)=xm
            zpym=zpm*ym
c           at(iz)=at(iz)+xm
c           zt(iz)=zt(iz)+zpym
c
c--abar, zbar, benuc, are the average baryon number, proton number 
c  number, and binding energy per nucleon (in Mev/nuc).
c
            ytst=ytst+ym
            atst=atst+xm
            ztst=ztst+zpym
            benuc=benuc+bim*ym
          ENDDO
c--normalize yeh
          if (xh.lt.1d-10) then
            xh=0.0d0
            yeh=0.0d0
          else
            yeh=yeh/xh 
          endif         
c         Do 41 i=1,32
c           fye=zt(i)/at(i)
c           write(15,2050) i,fye,zt(i),at(i)
c  41     Continue
          abar=atst/ytst
          zbar=ztst/ytst
c
c--convert from binding energy per nucleon to ergs per gram
c--with ye correction
          ubind=-9.616221376d17*(benuc+ye*0.783d0)
c         write(10,2000)(y(i),i=1,nmax) 
c         write(10,2100)atst,ztst
      Endif
      Return
  999 Format(a5,6(1x,1pd10.3))
 1000 FORMAT(3x,d9.2,6x,d9.2)
 1100 FORMAT(5x,d11.4,8x,d10.3,5x,i3,8x,i3)
 1200 FORMAT(4x,d13.6,5x,d13.6,4x,d11.4)
 1300 FORMAT(5x,d9.2,5x,i3,4x,f4.1,4x,f4.1,4x,f4.1,4x,f4.2,4x,f4.2,
     *4x,f4.2)
 2000 FORMAT (4(1pd15.7))
 2050 Format(i2,3(1x,1pd10.3))
 2100 FORMAT(1x,'Mass Conserv=',1pd15.7,' and Ye=',1pd15.7)
 2200 Format(1x,6d11.3)
 2300 Format(a5,1x,d9.2,1x,a5,d9.2,1x,a5,d9.2,1x,a5,1x,d9.2)
 3000 FORMAT(1x,'Convergence failure; Det=0 in step',I3,' at T9=',
     *1pd9.2,' and density=',1pd9.2)
 3100 FORMAT(1x,'Conv successful for T9=',1pd11.4,
     *' rho=',1pd11.4,' ye=',1pd11.4)
 3200 FORMAT(1x,'Convergence failure; does not converge in',I3,
     *' steps at T9=',1pd9.2,' and density=',1pd9.2)
 4000 FORMAT(1x,'k=',i3,' yp=',1pd16.8,' yn=',1pd16.8)
 4100 Format(1x,' delyp=',1pd13.6,' delyn=',1pd13.6,' det=',1pd13.6)
 4200 Format(1x,' dfp=',1pd11.4,' dfn=',1pd11.4,' dgp=',1pd11.4,' dgn=', 
     * 1pd11.4)     
 4250 Format(1x,' dfp2=',1pd11.4,' dfn2=',1pd11.4,' dgp2=',1pd11.4,
     * ' dgn2=',1pd11.4)     
 4300 format(1x,'f=',1pd11.4,'g=',1pd11.4)
 4400 Format(1x,f4.1,1x,f4.1,1x,f7.4,1x,1pd10.3,1x,'scr')
      End
c
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
  100 CONTINUE                                          
      t9i(24)=t9i(24)*10.                                                       
      DO 110 n=1,nyab
        read( 8,1000) inam(n) 
  110 CONTINUE                  
c--------------------------------------------------------------------------
c  This subroutine reads, from the file netwinv3, the nuclear data which will
c  be needed for later calcultions.  This data includes the atomic number,   
c  the number of protons and neutrons, and the binding energy (calculated)
c  from the tabulated mass excess.  Also the tabulation of the  partition
c  function, g, are read in for later interpolation. Once the set of nuclear
c  data is read in, it is assigned to the proper nuclei.
c-------------------------------------------------------------------------- 
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
  120 CONTINUE    
      bi(1)=0.d0
      bi(2)=0.d0                                                    
c     write(15,3000) (inam(i),za(i),zp(i),zn(i),bi(i),i=1,nyab)        
  999 Format(6(1x,1pd12.5))
 1000 format(a5)                                                               
 1010 format(24i3)                                                    
 2000 format(a5,f12.3,2i4,f6.1,f10.3)                                 
 2001 format(8f9.2)                                                   
 3000 format(1x,a5,3f15.2,f15.3)                                                
      RETURN                                                                    
      END                           
c
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
c--------------------------------------------------------------------------
c  The first step in calculating alph is to interpolate the correct g, 
c  the partition function, for the temperature.                           
c-------------------------------------------------------------------------- 
      do j=1,24
         i=j 
         if(t9i(j).gt.t9)go to 225
      enddo
  225 continue
c
      bkt=bok*t9 
      dbkt=1.d0/bkt
c     tmp=1.5d0*dlog(((avn*rho)**(2.d0/3.d0)*2.d0*pi*hbr**2)/(bkt*amu))
      tmp=dlog(rho*avn)+c1-1.5d0*dlog(T9)
c     If(t9.eq.5.) Then
c         write(15,*)rho,avn,hbr,pi,bkt,tmp
c     Endif
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
c                    
      SUBROUTINE nsescreen(T9,rho,ye,h,dhdp,dhdn,dhdt9,iscrn)
      implicit double precision (a-h,o-z)
      dimension h(37),dhdp(37),dhdn(37),dhdt9(37)
c
      third=1.d0/3.d0
c
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
C-----STRONG SCREENING BY ITOH ET AL.(1990)
        GNZ=ZZm*GNP
        dGNZdT9=ZZm*dGNPdT9
c       GNZ=GNP*ZZ*2.d0/(Z13+Z23)                                           
        EM2=A1*A2*2.d0/(A1+A2)                                              
        TAU=3.3722d0*(EM2*ZZ*ZZ/T9)**third
        dTAUdT9=TAU/(-3.d0*T9)
        GT=GNZ/TAU                                                        
        GT3=3.d0*GT
        GT33=GT3*GT3*GT3
        GT36=GT33*GT33
        GT312=GT36*GT36
        FNUM=.0455d0*GT3+.348d0*GT33+9.49d0*GT36-.123d0*GT312+
     #    .101d0*GT312*GT3
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
c           write(*,*) 'H',h0,dh0dt9
        Endif                                    
        FSR=GNZ*(1.25d0-1.71d0*GT)*0.333333333333333333d0
        FST=-GNZ*(1.25d0-1.425d0*GT)  
  30    Continue
c-----
c  Add succeeding screening factors
        Htot=Htot+h0
        dHtdp=dHtdp+dhdyp
        dHtdn=dHtdn+dhdyn
        dHtdT9=dHtdT9+dH0dT9
        h(j+1)=Htot
        dhdp(j+1)=dHtdp
        dhdn(j+1)=dHtdp
        dhdT9(j+1)=dHtdT9
        scrfct=exp(h(j+1))
c       If(T9.eq.Tfnl.and.rho.eq.rhofnl) Then
c           write(15,999) 'scrn',z1,z2,h0,dh0dp,dh0dn
c           tstscr=exp(h0)
c           write(15,*) 'TAU',tau,gnz                              
c           write(15,999) 'total',h(j+1),dhdp(j+1),dhdn(j+1)
c           write(14,999) 'scrn',z1,z2,h0,tstscr,scrfct
c       Endif
  100 Continue                        
      RETURN 
  999 Format(a5,6(1x,1pd10.3)) 
      END 
c
      subroutine newnorm(x,za,n1,n,zp,zn,ye)
      implicit double precision (a-h,o-z)
      dimension x(n),za(n),zp(n),zn(n)
      a=0.0d0
      b=0.0d0
      c=0.0d0
      d=0.0d0
c-----------------------------------------------------c
c   summation over loop n1 covers read-in abundances  c
c   if mass fractions ommit za(i)                     c  
      Do 10 i=1,n1
        a=a+zn(i)*x(i)
        b=b+zp(i)*x(i)
        c=c+zn(i)*zp(i)*x(i)/za(i)
        d=d+(zp(i)**2)*x(i)/za(i)
c       write(6,*) x(i),za(i),zp(i),zn(i),a,b
   10 Continue
c     write(15,*) 'newnorm',a,b,c,d
c-----------------------------------------------------c
      n2=n1+1
      s2=0.d0
      if(n2.gt.n) goto 18
c-----------------------------------------------------
c   takes solar abundances from n2 to n
c-----------------------------------------------------
      do 15 i=n2,n
        a=a+zn(i)*x(i)                 
        b=b+zp(i)*x(i)
        c=c+zn(i)*zp(i)*x(i)/za(i)
        d=d+zp(i)**2*x(i)/za(i)
   15 Continue
   18 beta=(ye*a-c)/(a*d-b*c)
      alph=(1-beta*b)/a
c-----------------------------------------------------c
      Do 20 i=1,n1
        x(i)=x(i)*(alph*zn(i)+beta*zp(i))/za(i)
   20 Continue
c-----------------------------------------------------c
      return
      end 


