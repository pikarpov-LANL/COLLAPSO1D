      program standalonemix
c
      implicit none
c
      double precision rho(4000),vr(4000),pr(4000),
     $     vsound(4000),menc(4000),r(4000),hp(4000)
      double precision bvf(4000),geff(4000),lam(4000)
      double precision vturb2(4000),lmix(4000)
      double precision alphalam,alphae,alphaye,alphak
      double precision dm(4000)
      double precision gg
c
c-- derivatives
      double precision dvturb(4000)
      double precision dvdr,drhodr,dprdr,dr
      double precision dpdt(4000),dudt(4000)
      double precision dt,dt0,time
      double precision ktilde(4000),dktilde(4000),
     $     a3(4000)
c
      integer ncell,i,j,iagain,nmax,nmaxp,iend,ijnk
      integer jd,idump
c
      character*20 outfile
c
      print *, 'dump number'
      jd=2000
      read *, idump
      print *, idump
      alphalam=0.1
      alphae=0.1
      alphaye=0.1
      alphak=0.1
c
      gg=6.67d-8
c
c      call readmesa(rho,vr,pr,vsound,menc,r,ncell)
      call readcollapse(rho,vr,pr,vsound,menc,r,hp,
     $     idump,ncell,outfile)
c     
      open(69,file=outfile)
      iend=0
      do i=1,ncell
         vturb2(i)=1.d2
         ktilde(i)=1e-2
      end do

      time=0.d0
      do j=1,jd
         dt=1.d1
         nmax=0
         do i=1,ncell
            if (i.eq.1) then
               dm(i)=menc(i)
            else 
               dm(i)=menc(i)-menc(i-1)
            end if
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
            geff(i)=gg*menc(i)/r(i)**2
c+vr(i)*dvdr
            bvf(i)=-geff(i)/rho(i)*(drhodr-
     $           dprdr/vsound(i)**2)
            bvf(i)=max(0.,bvf(i))
            lam(i)=alphalam*pr(i)/rho(i)/(gg*menc(i)/r(i)**2)
            lmix(i)=alphalam*hp(i)
c
            a3(i)=0.39*lam(i)*ktilde(i)**0.5*drhodr/rho(i)
c
c--calculating turbulent velocity evolution
c
            if (i.eq.1.or.geff(i).lt.0) then
               dvturb(i)=0.d0
               dktilde(i)=0.d0
            else
               dvturb(i)=r(i)**2*rho(i)*(vturb2(i)*vr(i)-
     $              alphak*lmix(i)*(vturb2(i)-vturb2(i-1))/dr)-
     $              r(i-1)**2*rho(i-1)*(vturb2(i-1)*vr(i-1)-
     $              alphak*lmix(i-1)*(vturb2(i)-vturb2(i-1))/dr)
               dvturb(i)=dvturb(i)/dr/r(i)**2
               dktilde(i)=a3(i)/rho(i)*dprdr*dm(i)
               if (i.eq.20) then
                  print *,dvturb(i) 
               end if


cdvturb(i)-
               dvturb(i)=-rho(i)*vturb2(i)*dvdr+
     $              rho(i)*dsqrt(vturb2(i))*bvf(i)*lmix(i)-
     $              rho(i)*dsqrt(vturb2(i))**3/lmix(i)
               if (i.gt.10) then
                  write(70,102)i,menc(i)/2.d33,dvturb(i)/vsound(i)**2.0,
c     $                 -rho(i)*vturb2(i)*dvdr/vsound(i)**2,
c     $              rho(i)*dsqrt(vturb2(i))*bvf(i)*lmix(i)/vsound(i)**2,
c     $              rho(i)*dsqrt(vturb2(i))**3/lmix(i)/vsound(i)**2,
     $                 dktilde(i)/vsound(i)**2/dm(i)
c                  write(70,102) i,bvf(i),lmix(i),rho(i),vturb2(i)
               end if
               
            end if

            dvturb(i)=max(0.d0,dvturb(i))
            if (vturb2(i).gt.0.95*vsound(i)**2) then
               dvturb(i)=0.
               nmax=nmax+1
            end if
            dt0=0.1*(vsound(i)**2/(abs(dvturb(i))+1.d-10))
            dt0=min(dt0,(vturb2(i)/(abs(dvturb(i))+1.d-10)))
            if (dt0.lt.dt) then
               dt=dt0
            end if
            if (i.eq.20) then
               print *, 'final',vturb2(i),dvturb(i),vsound(i)**2,time
               print *, 'livescu',ktilde(i),dktilde(i),a3(i)
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
         time=time+dt
         do i=1,ncell
            vturb2(i)=vturb2(i)+dvturb(i)*dt
            ktilde(i)=ktilde(i)+dktilde(i)*dt
            ktilde(i)=max(1e-8,ktilde(i))

            if (j.eq.jd.or.time.gt.0.001) then
               if (i.eq.1) write(69,102) j, time, dt
               write(69,102) i,menc(i)/2.d33,dsqrt(vturb2(i))/vsound(i),
     $              dvturb(i)/vsound(i)**2,
     $              ktilde(i),
     $              dktilde(i)/vsound(i)**2/dm(i),
     $              ktilde(i)/vsound(i)**2/dm(i)
               iend=1
            end if
         end do
         if (iend.eq.1) stop
      end do
 102  format(I4,9(1pe12.4))
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


      subroutine readcollapse(rho,vr,pr,vsound,menc,r,hp,
     $     idump,ncell,outfile)
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
     $     vsound(4000),menc(4000),r(4000),temp(4000),
     $     ye(4000),hp(4000)
      character*30 infile
      character*1 cjnk
      character*20 outfile
      double precision djnk,tcol
      integer i,ijnk,ncell,izone,j,idump
c
      do i=1,idump
         open(42,file='input')
         read(42,*) ncell,infile,outfile
      end do
      ncell=ncell-1
      close(42)
      open(42,file=infile)
      read(42,*) tcol
      do i=1,ncell
         read(42,*)izone,menc(i),r(i),rho(i),vr(i),ye(i),
     $        pr(i),vsound(i)
         menc(i)=2.d33*menc(i)
      end do
      do i=1,ncell-50
         do j=i+1,ncell
            hp(i)=r(j)-r(i)
            if (pr(j).lt.pr(i)/2.718218) goto 20
         end do
 20      continue
      end do
      do i=ncell-51,ncell
         hp(i)=hp(ncell-50)
      end do
      close(42)
      return
      end

