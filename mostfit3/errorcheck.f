      program errorchecker
      implicit none
      integer nmax,nunit,npt,nbin,i,j,k,bin,nlpt,ii,jj
      parameter(nmax=500000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),skystd(nmax),
     .   skymin,skymax,binsize,bdatax(nmax),bdatay(nmax),lsky(nmax),
     .   lerr(nmax),std,stdev,mean,lmag(nmax),meansky,phi,phase(nmax),
     .   lphase(nmax),Per,dph,Pi,zphase,Eerr,rmavg,errfs(nmax),
     .   errfp(nmax),btemp(nmax)
   
      double precision dtime(nmax),ldtime(nmax),ktime(nmax)
      character*80 filename
      Per=3.52474859
      Pi=acos(-1.0)
      
      dph=0.02
      phi=-0.6820788 !phase offset in radians
      zphase=0.5-phi/(2.0*Pi)
      
      write(6,*) "Enter filename"
      read(5,500) filename
 500  format(A80)
 
      nunit=10
      open(unit=nunit,file=filename,status='old',err=901)
      
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)     
      close(10)
      
      call phasept(npt,time,mag,phase,Per)
      
      j=0
      do 5 i=1,npt
         if((sky(i).gt.0.0).and.(mfield(i).gt.20000.0))then
c            if((mag(i).gt.-0.1).and.(mag(i).lt.0.1))then
c               if((phase(i).lt.zphase-dph).or.
c     .          (phase(i).gt.zphase+dph))then
                  j=j+1
c                  if((sky(i).gt.120.0).and.(sky(i).lt.600.0))
c     .            merr(i)=merr(i)*(sky(i)/1000.0+0.89) 
c                  if(sky(i).ge.600.0) merr(i)=merr(i)*1.5
                  lsky(j)=log(sky(i))
                  lerr(j)=merr(i)
                  lmag(j)=mag(i)
                  lphase(j)=phase(i)
                  ldtime(j)=dtime(i)
c               endif
c            endif
         endif
 5    continue
      nlpt=j
      
      call avgrm(nlpt,lmag,rmavg)
      call findrange(nlpt,lsky,skymin,skymax)
      write(6,*) "sky range: ",skymin,skymax
      
      nbin=100
C     get size of bin.
      binsize=(skymax-skymin)/real(nbin-1)
      
      call pgopen('?')
      ii=0
      do 10 i=1,nbin
         k=0
         meansky=0.
         do 11 j=1,nlpt
            bin=int((lsky(j)-skymin)/binsize)+1
            if(bin.eq.i)then
               k=k+1
               bdatax(k)=lmag(j)
               bdatay(k)=lerr(j)
               ktime(k)=ldtime(j)
               meansky=meansky+exp(lsky(j))
            endif
 11      continue
         Eerr=0.
         do 12 j=1,k
c            Eerr=Eerr+1.0/bdatay(j)
            Eerr=Eerr+bdatay(j)*bdatay(j)
 12      continue
c         Eerr=sqrt(real(k))/Eerr
         Eerr=Sqrt(Eerr)/Sqrt(real(k))
         if(k.gt.100) then
            ii=ii+1
c            mean=0.0

            do 14 jj=1,k
               btemp(jj)=bdatax(jj)
 14         continue
            jj=k
            call sigcut(jj,btemp,3.0)
            std=stdev(jj,btemp,mean)
            
            errfs(ii)=meansky/real(k)
            errfp(ii)=std/Eerr
            write(6,*) errfs(ii),errfp(ii),std,Eerr
            call plotbins(k,bdatax,mean,std,Eerr)
c            do 13 j=1,k
c               if(abs(bdatax(j)).gt.4.0*std)then
c                  write(6,501) ktime(j),bdatax(j)
c 501              format(F13.8,1X,2(F9.6,1X))
c               endif
c 13         continue
         endif
 10   continue
      call pgclos()

      call errorcor(npt,merr,sky,ii,errfs,errfp)
      rmavg=0.
      filename="errcor.dat"
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)

      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorcor(npt,merr,sky,ii,errfs,errfp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,ii,nmax,i
      parameter(nmax=500 000)
      real merr(npt),sky(npt),errfs(ii),errfp(ii),yp1,yp2,y2(nmax),xx,
     .   yy
      
      yp1=1.0e30
      yp2=1.0e30
      call spline(errfs,errfp,ii,yp1,yp2,y2)
      
      do 10 i=1,npt
         if(sky(i).gt.140.0)then
            xx=sky(i)
            call splint(errfs,errfp,y2,ii,xx,yy)
            merr(i)=merr(i)*yy
         endif
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotbins(npt,xdata,mean,std,Eerr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,nmax,nplot,i
      parameter(nmax=500 000,nplot=1000)
      real xdata(npt),std,bdatax(nmax),bdatay(nmax),bmax,xb1,xb2,yb1,
     .   yb2,fc,datamin,datamax,mean,plotx(nplot),ploty(nplot),rgauss,
     .   area,Eerr
      
      nbin=real(npt)/20.0
      datamin=mean-4.0*std
      datamax=mean+4.0*std
c      call findrange(npt,xdata,datamin,datamax)
      call bindata(nbin,npt,xdata,bdatax,bdatay,datamin,datamax,
     .      bmax)  !bin the data.

      call pgpage() !prepare for multiple pages
      call windowsetup(xb1,xb2,yb1,yb2) !setup a square plot
      call pgvport(xb1,xb2,yb1,yb2) !make a square plotting surface
      fc=2.*(datamax-datamin)/real(nbin)
c      fc=0.
c         write(6,*) "fc:",fc
      call pgwindow(datamin+fc,datamax+fc,0.0,bmax+0.1*bmax)!set window scale
      call pgslw(5)
      call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)!boxes and ticks
      call pglabel("Measurements (mag)","# of Values","") !label axes     
      call pgbin(nbin,bdatax,bdatay,.false.)
      
      area=0.
      do 11 i=1,nbin
         area=area+bdatay(i)*fc/2.0
 11   continue
      
      do 10 i=1,nplot
         plotx(i)=(datamax-datamin)*real(i)/real(nplot)+datamin
         ploty(i)=area*rgauss(std,mean,plotx(i))
 10   continue
 
      call pgsci(2)
      call pgline(nplot,plotx,ploty)
      call pgsci(1)
      
      do 12 i=1,nplot
         plotx(i)=(datamax-datamin)*real(i)/real(nplot)+datamin
         ploty(i)=area*rgauss(Eerr,mean,plotx(i))
 12   continue
 
      call pgsci(3)
      call pgline(nplot,plotx,ploty)
      call pgsci(1)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function rgauss(sig,mu,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real sig,mu,x,Pi,temp(2)
      Pi=3.14159265
      
      temp(1)=sig*sqrt(2.0*Pi)
c      temp(1)=1.0
      temp(2)=(x-mu)*(x-mu)/(2*sig*sig)
      rgauss=exp(-temp(2))/temp(1)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine avgrm(npt,mag,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Removes average value from data points
C     This can help reduce numerical overflow for large data sets 
      implicit none
      integer npt
      real mag(npt),rmavg

      integer i
      real ave,adev,sdev,sigma2,skew,curt

      ave=0.0
      do 9 i=1,npt
         ave=ave+mag(i)
 9    continue
      ave=ave/real(npt)

      do 10 i=1,npt
         mag(i)=mag(i)-ave
 10   continue

      rmavg=ave
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,mag,phase,period)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      real time(npt),mag(npt),phase(npt),period

      integer i
      real temp

      do 10 i=1,npt
c         temp=time(i)/period
c         phase(i)=temp-int(temp)
         temp=time(i)!-time(1)
         phase(i)=temp/period-int(temp/period)
c         write(6,*) time(i),phase(i)

c         phase(i)=phase(i)+0.5
c         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigcut(npt,pts,sig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i,j,k,npt1,niter
      parameter(nmax=500000)
      real pts(npt),tmp1(nmax),tmp2(nmax),std,mean,stdev,sig
      
      do 10 i=1,npt
         tmp1(i)=pts(i)
 10   continue      
      npt1=npt
      
      niter=1
      
      do 22 k=1,niter
         std=stdev(npt1,tmp1,mean)
         j=0
         do 20 i=1,npt
            if(abs(tmp1(i)-mean).lt.sig*std)then
               j=j+1
               tmp2(j)=tmp1(i)
            endif
 20      continue
         do 21 i=1,j
            tmp1(i)=tmp2(i)
 21      continue
         npt1=j
 22   continue
 
      do 30 i=1,npt1
         pts(i)=tmp1(i)
 30   continue
      npt=npt1
      
      return
      end
 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason F. Rowe (2005)
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i,ntmp
      real pts(npt),mean,s,ep,p,var,sdev

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      ntmp=0
      do 10 i=1,npt
         if(pts(i).ne.99.9) then
            ntmp=ntmp+1
            s=pts(i)-mean
            ep=ep+s
            p=s*s
            var=var+p
         endif
 10   continue
      var=(var-ep**2/ntmp)/(ntmp-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(npt,pdata,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real pdata(npt),datamin,datamax
      
C     initialize datamin and datamax
      datamin=pdata(1)
      datamax=pdata(1)
      
      do 10 i=2,npt
         datamin=min(pdata(i),datamin)
         datamax=max(pdata(i),datamax)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i)*binsize
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowsetup(xb1,xb2,yb1,yb2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      setup plotting area
      real xb1,xb2,yb1,yb2,xs1,xs2,
     .   ys1,ys2,pratio
      
C     set text size      
      call pgsch(1.0)
C     inquire window size
      call pgqvsz(2,xs1,xs2,ys1,ys2)
      pratio=(xs2-xs1)/(ys2-ys1)
C     if we have a landscape surface
      if(xs2-xs1.gt.ys2-ys1)then
C        start 10% from the edge
         yb1=(ys2-ys1)*0.10+ys1
C        use 60% of surface available
         yb2=(ys2-ys1)*0.70+yb1
C        now normalize to 1
         yb1=yb1/(ys2-ys1)         
         yb2=yb2/(ys2-ys1)
C        now use same size in other dimension      
         xb1=yb1
         xb2=yb1+(yb2-yb1)/pratio
      else
C        start 10% from the edge
         xb1=(xs2-xs1)*0.10+xs1
C        use 60% of surface available
         xb2=(xs2-xs1)*0.60+xb1
C        now normalize to 1
         xb1=xb1/(xs2-xs1)         
         xb2=xb2/(xs2-xs1)
C        now use same size in other dimension      
         yb1=xb1
         yb2=xb1+(xb2-xb1)*pratio
      endif
      
      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE spline(x,y,n,yp1,ypn,y2) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n,NMAX 
      REAL yp1,ypn,x(n),y(n),y2(n) 
      PARAMETER (NMAX=400000) 
      INTEGER i,k 
      REAL p,qn,sig,un,u(NMAX) 
      if (yp1.gt..99e30) then 
         y2(1)=0. 
         u(1)=0. 
      else 
         y2(1)=-0.5 
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 
      do 11 i=2,n-1 
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
         p=sig*y2(i-1)+2. 
         y2(i)=(sig-1.)/p 
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     .        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
 11   continue
      if (ypn.gt..99e30) then 
         qn=0. 
         un=0. 
      else 
         qn=0.5 
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.) 
      do 12 k=n-1,1,-1 
         y2(k)=y2(k)*y2(k+1)+u(k) 
 12   continue
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE splint(xa,ya,y2a,n,x,y) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n 
      REAL x,y,xa(n),y2a(n),ya(n) 
      INTEGER k,khi,klo 
      REAL a,b,h 
      klo=1 
      khi=n 
 1    if (khi-klo.gt.1) then 
         k=(khi+klo)/2 
         if(xa(k).gt.x)then 
            khi=k 
         else 
            klo=k 
         endif 
      goto 1 
      endif 
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) pause 'bad xa input in splint'  
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi)+ 
     .     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6. 
      return 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gasdev(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
 1       v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END