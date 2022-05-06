      program standalonemix
c
      implicit none
c
      double precision rho(4000),vr(4000),pr(4000),
     $     vsound(4000),menc(4000),r(4000)
      double precision bvf(4000),geff(4000),lam(4000)
      double precision vturb2(4000),lmix(4000)
      double precision alphalam,alphae,alphaye,alphak
      double precision gg
c
c-- derivatives
      double precision dvturb(4000)
      double precision dvdr,drhodr,dprdr,dr
      double precision dpdt(4000),dudt(4000)
      double precision dt,dt0,time
c
      integer ncell,i,j,iagain,nmax,nmaxp
c
      alphalam=0.1
      alphae=0.1
      alphaye=0.1
      alphak=0.1
c
      gg=6.67d-8
c
      call readmesa(rho,vr,pr,vsound,menc,r,ncell)
c     
      do i=1,ncell
         vturb2(i)=1.d0
      end do

      time=0.d0
      do j=1,5000
         dt=1.d10
         nmax=0
         do i=1,ncell
            if (i.eq.1) then
               dr=r(i)
               dvdr=vr(i)/dr
               drhodr=0.d0
               dprdr=0.d0
            else
               dr=(r(i)-r(i-1))
               dvdr=(vr(i)-vr(i-1))/dr
               drhodr=(rho(i)-rho(i-1))/dr
               dprdr=(pr(i)-pr(i-1))/dr
            end if
            geff(i)=gg*menc(i)/r(i)**2+vr(i)*dvdr
            bvf(i)=-geff(i)/rho(i)*(drhodr-
     $           dprdr/vsound(i)**2)
            lam(i)=alphalam*pr(i)/rho(i)/geff(i)
            lmix(i)=alphalam*pr(i)/rho(i)/geff(i)
c
c--calculating turbulent velocity evolution
c
            if (i.eq.1) then
               dvturb(i)=0.d0
            else
               dvturb(i)=r(i)**2*rho(i)*(vturb2(i)*vr(i)-
     $              alphak*lmix(i)*(vturb2(i)-vturb2(i-1))/dr)-
     $              r(i-1)**2*rho(i-1)*(vturb2(i-1)*vr(i-1)-
     $              alphak*lmix(i-1)*(vturb2(i)-vturb2(i-1))/dr)
               dvturb(i)=dvturb(i)/dr/r(i)**2
               dvturb(i)=dvturb(i)-rho(i)*vturb2(i)*dvdr+
     $              rho(i)*dsqrt(vturb2(i))*bvf(i)*lmix(i)-
     $              rho(i)*dsqrt(vturb2(i))**3/lmix(i)
            end if
            dvturb(i)=max(0.d0,dvturb(i))
            if (vturb2(i).gt.0.95*vsound(i)**2) then
               dvturb(i)=0.
               nmax=nmax+1
            end if
            dt0=0.1*(vsound(i)**2/(dvturb(i)+1.d-10))
            if (dt0.lt.dt) then
               dt=dt0
            end if
            if (dt.lt.0) then
               print *, i,dt,vsound(i),dvturb(i),(dvturb(i)+1.d-10)
               read *, iagain
            end if
c
c- momentum equation
c
            if (i.eq.1) then
               dpdt(i)=0.d0
            else
               dpdt(i)=rho(i)*vturb2(i)
            end if
c
c--energy equation
c
            if (i.eq.1) then
               dudt(i)=0.
            else
               dudt(i)=rho(i)*vturb2(i)**1.5/lmix(i)
            end if
         end do
         if (nmax.gt.nmaxp)  print *, j,time/3.15d7,nmax
         nmaxp=nmax
         do i=1,ncell
            vturb2(i)=vturb2(i)+dvturb(i)*dt
            time=time+dt
            if (j.eq.1000) then
               write(69,102) j,i,vturb2(i)/vsound(i)**2,
     $              dpdt(i)/pr(i)
            end if
         end do
      end do
 102  format(2(I5),2(1pe10.2))
      end

      subroutine readmesa(rho,vr,pr,vsound,menc,r,ncell)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                               c
c This subroutine reads in the stellar data from Mesa           c
c                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      double precision rho(4000),vr(4000),pr(4000),
     $     vsound(4000),menc(4000),r(4000),temp(4000)
      character*30 infile
      character*1 cjnk
      double precision djnk
      integer i,ijnk,ncell
c
      open(42,file='inputstar')
      read(42,*) infile
      close(42)
      open(42,file=infile)
      read(42,*) ijnk
      read(42,*) cjnk
      read(42,*) ijnk,ncell
      read(42,*) ijnk
      read(42,*) cjnk
      do i=ncell,1,-1
         read(42,*)ijnk,temp(i),rho(i),pr(i),r(i),
     $        djnk,djnk,vr(i),djnk,djnk,
     $        vsound(i),djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        djnk,djnk,djnk,djnk,djnk,
     $        menc(i),djnk,djnk,djnk,djnk
         temp(i)=10**temp(i)
         rho(i)=10**rho(i)
         pr(i)=10**pr(i)
         r(i)=7.d10*10**r(i)
         menc(i)=2.d33*menc(i)
      end do
      close(42)
      return
      end

