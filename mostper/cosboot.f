      program cosinefitbootstrap
      implicit none
      integer npt,nmax,filet,ma,nfmax,nstar,i,nfit,now(3),seed,niter,j,
     .	iargc
      parameter(nmax=500000,nfmax=2000,niter=10000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),ftotal(nmax),
     .   aflux(nmax),nap(nmax),etime(nmax),skystd(nmax),mintime,
     .   mfield(nmax),a(nfmax),rmavg,aerr(nfmax),peak(nmax),sig,
     .   x(nmax),y(nmax),z(nmax),dumr,ran2,ao(nfmax)
      double precision dtime(nmax)
      character*80 filename
      common /data/ x,y,z,npt,nfit

C     Random number creation (so we get a different seed each time)
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)

      if(iargc().lt.3) goto 903
      
      
      call getarg(1,filename)
c      write(6,*) "Enter MOST photometry data file"
c      read(5,*) filename
      open(unit=10,file=filename,status='old',err=901)

      call getarg(2,filename)
c      write(6,*) "Enter Frequency database"
c      read(5,*) filename
      open(unit=11,file=filename,status='old',err=901)     
      
      call getarg(3,filename)
c      write(6,*) "Enter name for bootstrap results"
c      read(5,*) filename
      open(unit=12,file=filename)
      
      filet=0
      call readdata(npt,time,dtime,mag,merr,sky,xc,yc,fx,fy,fxy,tboard,
     .   mfield,ftotal,aflux,nap,etime,skystd,filet)
c      write(6,*) "Number of datapoints :", npt
      
      call avgrm(npt,mag,rmavg)
      call avgrm(npt,mag,rmavg)
      
      call readfreqs(ma,nstar,ao,mintime)
      
      do 10 i=1,npt
         time(i)=time(i)-mintime
 10   continue
 
c      do 5 i=1,npt
c         write(6,*) time(i),mag(i),merr(i)
c 5    continue 
 
      
      sig=3.0
      do 11 i=1,1
         call sigclip(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,fxy,
     .     tboard,peak,ftotal,aflux,nap,etime,rmavg)
 11   continue
      
      do 12 i=1,niter
         write(6,*) "Iteration #:",i
         do 13 j=1,ma
            a(j)=ao(j)  !original solution
 13      continue
         call bootdata(seed,npt,time,mag,merr,x,y,z)
         call cosinefit(a,rmavg,ma,nstar,aerr)
         call freqout(ma,nstar,a)
 12   continue
 
      close(10)
      close(11)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 903  write(6,*) "Usage: cosboot <photometry> <frequencies> <bootstrap>"
      write(6,*) "<photometry>: name of MOST photometry file"
      write(6,*) "<frequency>: frequency solution file"
      write(6,*) "<bootstrap>: filename for bootstrap output"
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bootdata(seed,npt,time,mag,merr,x,y,sig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer seed,npt,i,j
      real time(npt),mag(npt),merr(npt),x(npt),y(npt),sig(npt),ran2

C     resample data
      do 10 i=1,npt
         j=ran2(seed)*npt+1.0
         x(i)=time(j)
         y(i)=mag(j)
         sig(i)=merr(j)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine freqout(ma,nstar,a)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ma,nstar,nfmax,i,j,k,cc,b,nstarmax
      parameter(nfmax=2000,nstarmax=100)
      real a(nfmax),twopi,pi,aerr(nfmax),sns(nstarmax)

      pi=3.141592654
      twopi=2.0*pi

      j=(ma-nstar)/nstar
      do 38 k=1,nstar
         cc=(j+1)*(k-1)+2
         a(cc-1)=a(cc-1)/twopi
         do 39 i=cc,cc+j-2,2
c            if(a(i).lt.0.0) then
c               a(i)=abs(a(i))
c               a(i+1)=a(i+1)+pi
c            endif
            b = int(a(i+1)/twopi)
            a(i+1)=a(i+1)-b*twopi
c            if (a(i+1).lt.0) a(i+1)=twopi+a(i+1)
 39      continue
 38   continue

      write(12,500) (a(i),i=1,ma)
      
c      j=(ma-nstar)/nstar
c      do 10 k=1,nstar
c         cc=(j+1)*(k-1)+2
c         write(12,500) a(cc-1)/twopi,(a(i),i=cc,cc+j-2,2),
c     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
c 10   continue

 500  format(300(E14.7,1X))
c 500  format(20(F13.6,1X))
 501  format(A3,1X,F11.4)
 502  format(I4,1X,I4)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      real pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigclip(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,
     .     fxy,tboard,peak,ftotal,aflux,nap,itime,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=500000,nitmax=50)
      real time(npt),mag(npt),merr(npt),sig,std,stdev,tmp1(nmax),
     .     tmp2(nmax),tmp3(nmax),mean,dum,tmp4(nmax),tmp5(nmax),
     .     tmp6(nmax),tmp7(nmax),tmp8(nmax),tmp9(nmax),tmp10(nmax),
     .     tmp11(nmax),tmp12(nmax),sky(npt),xc(npt),omean,
     .     yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),peak(npt),rmavg,
     .     tmp13(nmax),tmp14(nmax),tmp15(nmax),tmp16(nmax),ftotal(nmax),
     .     aflux(nmax),nap(nmax),itime(nmax)


C     watch out for infinite loops
      niter=0

C     find mean of data set
      ntmp=0
      mean=0.
      do 4 i=1,npt
         if(mag(i).lt.90.0) then
            ntmp=ntmp+1
            mean=mean+mag(i)
         endif
 4    continue
      write(6,*) mean,ntmp
      mean=mean/real(ntmp)
C     find standard dev. of data set.
      std=stdev(npt,mag,mean)

      write(6,*) "sigclip:",npt,mean,std

C     count number of clipped points
 6    nclip=0
C     count number of new points
      ntmp=0
      do 10 i=1,npt
         dum=abs(mag(i)-mean)
         if(dum.lt.sig*std) then
            ntmp=ntmp+1
            tmp1(ntmp)=time(i)
            tmp2(ntmp)=mag(i)
            tmp3(ntmp)=merr(i)
            tmp5(ntmp)=sky(i)
            tmp6(ntmp)=xc(i)
            tmp7(ntmp)=yc(i)
            tmp8(ntmp)=fx(i)
            tmp9(ntmp)=fy(i)
            tmp10(ntmp)=fxy(i)
            tmp11(ntmp)=tboard(i)
            tmp12(ntmp)=peak(i)
            tmp13(ntmp)=ftotal(i)
            tmp14(ntmp)=aflux(i)
            tmp15(ntmp)=nap(i)
            tmp16(ntmp)=itime(i)
         else
            nclip=nclip+1
         endif
 10   continue
C     if nothing is clipped, no point in continuing
      if(nclip.eq.0) goto 15
      
C     save copy of old mean
      omean=mean

C     find mean of new data set
      mean=0.
      do 5 i=1,ntmp
         mean=mean+tmp2(i)
 5    continue
      mean=mean/real(ntmp)
C     if mean doesn't change, we're done.
      if(mean.eq.omean) goto 15

C     if mean doesn't change, we're done
      if(nclip.gt.0) goto 15 
C     find st. dev. of new data set
      std=stdev(ntmp,tmp2,mean)
C     now restart clipping
      niter=niter+1
C     if we loop too much.. get out. (probably an infinite loop)
      if(niter.gt.nitmax) goto 15
      goto 6

 15   do 20 i=1,ntmp
         time(i)=tmp1(i)
         mag(i)=tmp2(i)-mean
         merr(i)=tmp3(i)
         sky(i)=tmp5(i)
         xc(i)=tmp6(i)
         yc(i)=tmp7(i)
         fx(i)=tmp8(i)
         fy(i)=tmp9(i)
         fxy(i)=tmp10(i)
         tboard(i)=tmp11(i)
         peak(i)=tmp12(i)
         ftotal(i)=tmp13(i)
         aflux(i)=tmp14(i)
         nap(i)=tmp15(i)
         itime(i)=tmp16(i)
 20   continue
      npt=ntmp
C     update zero point
      rmavg=rmavg+mean

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(ma,nstar,a,mintime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ma,nstar,j,i,cc,k,nfmax
      parameter(nfmax=2000)
      real a(nfmax),mintime,pi,twopi,dumr
      
      do 12 i=1,nfmax
         a(i)=0.0
 12   continue
      
      pi=3.141592654
      twopi=2.0*pi
      
      read(11,*) mintime
      
      read(11,502) ma,nstar
      
      j=(ma-nstar)/nstar
      do 10 k=1,nstar
         cc=(j+1)*(k-1)+2
         read(11,500) a(cc-1),(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2)
         a(cc-1)=a(cc-1)*twopi
 10   continue
      
 500  format(20(E14.7,1X))      
 502  format(I4,1X,I4)

      return
      end
      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(npt,time,dtime,mag,merr,
     .   sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,
     .   skystd,filet)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     id - id of star (returned)
C     xcoo,ycoo - co-ordinates of star (returned)
C     time,mag,merr - the data and associated error.
      implicit none
      integer nmax
      parameter(nmax=500000)
      integer npt,i,filet,rn
      real time(nmax),mag(nmax),merr(nmax),c,amag(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),skystd(nmax),ntest
      double precision dtime(nmax)
   
C      read(10,*,end=30) id,xcoo,ycoo

      i=1

      if(filet.eq.0) then
   9     continue
  10     read(10,*,err=9,end=20) dtime(i),mag(i),merr(i),sky(i),
     .      xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .      ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
            time(i)=dtime(i)
c            merr(i)=0.001
            if(mfield(i).lt.20000.0) goto 10
c            if(mag(i).gt.0.05) goto 10
c            if(mag(i).lt.-0.05) goto 10
c            if((sky(i).lt.0.0).or.(sky(i).gt.2000.0)) goto 10
c            merr(i)=abs(merr(i))
            if(merr(i).le.0.0) goto 10
c           if(sky(i).gt.1500.0) goto 10
c         if((time(i).gt.1690.4).and.(time(i).lt.1690.6)) goto 10
c         if((time(i).gt.1693.8).and.(time(i).lt.1694.2)) goto 10
c         if((time(i).gt.1697.4).and.(time(i).lt.1697.7)) goto 10
c         if((time(i).gt.1700.9).and.(time(i).lt.1701.2)) goto 10
C         rn=ntest(mag(i))
C         if(rn.eq.1) goto 10
C         rn=ntest(sky(i))
C         if(rn.eq.1) goto 10
            i=i+1
         goto 10
 20      continue
      elseif(filet.eq.1) then
 29      continue
 30      read(10,*,err=29,end=40) dtime(i),mag(i),merr(i),sky(i)
C     .      xc(i),yc(i)!,fx(i),fy(i)
           time(i)=dtime(i)
           if(sky(i).gt.2000.0) goto 30
           if(mag(i).gt.80.0) goto 30
           merr(i)=abs(merr(i))
           if(merr(i).le.0.0) goto 30
           i=i+1
         goto 30
 40      continue
      endif

 505  format(F10.5,1X,3(F7.4,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F8.1)
C     11.6/13.8
 510  format(F13.8,1X,2(F9.6,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F9.2),1X,F8.3,1X,F6.2,1X,F8.4)
 511  format(F11.6,1X,2(F9.6,1X),F8.2,5(1X,F9.3))

      npt=i-1

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function ntest(rr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine is used to test for NANs and other nasty numbers.
      implicit none
      real rr,thi,tlow

C     this defines the numerical boundary of anything that might be 
C     exciting.

      thi=  99.9e30
      tlow=-99.9e30

      if((rr.gt.tlow).and.(rr.lt.thi))then
         ntest=0
      else
         ntest=1
      endif

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

      call moment(mag,npt,ave,adev,sdev,sigma2,skew,curt)

      do 10 i=1,npt
         mag(i)=mag(i)-ave
 10   continue

      rmavg=ave
      
      return
      end

ccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cosinefit(a,rmavg,ma,nstar,aerr)
CCCCccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     The one and only Fourier Decomposition routine.
C     per - fixed period for data
C     id,xcoo,ycoo - id and co-ordinates of object.
C     data is passed by common block /data/ 
C       x - time of observations
C       y - the observations
C       sig - the associated error
C       n - number of data points
      implicit none
      integer nmax,nfmax,nfit,nfitm,nstar,nstarmax,plot,endfit
      parameter(nmax=500000,nfmax=2000,nfitm=2000,nstarmax=100)
      real w(nstarmax), x(nmax),y(nmax),sig(nmax),a(nfmax),chisq,tx,ty,
     .     alpha(nfmax,nfmax),covar(nfmax,nfmax),alamda,avg,
     .     avgc,per,ty2(nmax),tx2(nmax),max,min,amp,
     .     pi,maxamp,rmavg,pers(nstarmax),arg,
     .     twopi,ochi,alamhi,chitol,aerr(nfmax)
      integer i,n,ia(nfmax),ma,nca, j, k,cc,j1,
     .     i2,niter
      character*50 filename
      logical loop
      common /data/ x,y,sig,n,nfit
c      external fit1a

      alamhi=100000.
      chitol=10e-6
      avg=0.

      maxamp=2.0
      pi=3.141592654
      twopi=2*pi

      loop=.true.
      
      call avgrm(n,y,rmavg)
      
      do 10 i=1,ma
         ia(i)=1.0
 10   continue
      nca=nfmax

      j=(ma-nstar)/nstar
      do 32 k=1,nstar
         cc=(j+1)*(k-1)+2
C     UNCOMMENT HERE FOR FIXED PERIODS
c         ia(cc-1)=0
 32   continue

      alamda=-1.0

      endfit=0.
      i2=0
      niter=275
c      do 20 i2=1,niter
      do while (loop)
         i2=i2+1
         ochi=chisq
         call mrqmin(x,y,sig,n,a,ia,ma,covar,alpha,nca,chisq,alamda,w,
     .        nstar)
c         if(i2.ge.niter) then
         if((i2.gt.niter).or.(endfit.gt.50.0).or.(alamda.gt.alamhi))then
c            write(6,*) "i2,endfit,alamda,niter,alamhi"
c            write(6,*) "converge?",i2,endfit,alamda,niter,alamhi
            alamda=0.0
            call mrqmin(x,y,sig,n,a,ia,ma,covar,alpha,nca,chisq,alamda,
     .           w,nstar)
            loop=.false.
         endif

C     correct average
         avgc=0.
         do 35 j1=1,n
            ty=0.
            j=(ma-nstar)/nstar
            do 36  k=1,nstar
               cc=(j+1)*(k-1)+2
               do 37 i=cc,cc+j-2,2
                  arg=real((i-cc+2)/2)*a(cc-1)*x(j1)+a(i+1)
                  ty=ty+a(i)*cos(arg)
 37            continue
 36         continue
            avgc=avgc+y(j1)-ty
 35      continue
         avgc=avgc/n
         do 16 j=1,n
            y(j)=y(j)-avgc
 16      continue
         avg=avg+avgc
         if(chisq.gt.ochi) then
            endfit=0
         elseif(abs(chisq-ochi).lt.chitol) then
            endfit=endfit+1
         endif
      enddo
 50   continue
      niter=i2

      amp=0.

      do 21 i=1,n
         y(i)=y(i)+avg
 21   continue

      j=(ma-nstar)/nstar
      do 23 k=1,nstar
         cc=(j+1)*(k-1)+2
         a(cc-1)=abs(a(cc-1))
         pers(k)=twopi/a(cc-1)
         w(k)=1.0/pers(k)
 23   continue

      max=-99.9e10
      min= 99.9e10

      do 22 i2=1,1000
         j=(ma-nstar)/nstar
         do 40 k=nstar,nstar
            cc=(j+1)*(k-1)+2
            per=twopi/a(cc-1)
            tx=real(i2-1)*per/1000.0+x(1)
            tx2(i2)=tx/per-int(tx/per)
            ty2(i2)=0.
            do 41 i=cc,cc+j-2,2
               arg=real((i-cc+2)/2)*a(cc-1)*tx+a(i+1)
               ty2(i2)=ty2(i2)+a(i)*cos(arg)
 41         continue
 40      continue
         ty2(i2)=ty2(i2)+avg
         if(ty2(i2).gt.max) max=ty2(i2)
         if(ty2(i2).lt.min) min=ty2(i2)
 22   continue

      amp = max-min
c      if(1.0/per.gt.5.0) then
c      write(6,*) "freq,btheta,avg,amp,rchi,niter,alamda"
c      write(6,502)1.0/per,btheta,avg+rmavg,amp,chisq/(n-1),niter,alamda
 502  format(F9.4,1X,F10.3,1X,F9.4,1X,F9.4,1X,F9.4,1X,I4,F9.4)
      
      j=(ma-nstar)/nstar
      do 42 k=1,nstar
         cc=(j+1)*(k-1)+2
         aerr(cc-1)=sqrt(covar(cc-1,cc-1))/twopi
         if(cc+j-2.eq.nfmax)write(6,*) "Increase nfmax to: ",cc+j-2
c         write(6,501) "Fre", a(cc-1)/twopi,aerr(cc-1)
         do 43 i=cc,cc+j-2,2
            aerr(i)=sqrt(covar(i,i))
            aerr(i+1)=sqrt(covar(i+1,i+1))
 43      continue
c         write(6,500) "Amp", (a(i)  ,i=cc,cc+j-2,2)
c         write(6,500) "err", (sqrt(abs(covar(i,i))),    i=cc,cc+j-2,2)
c         write(6,500) "Psi", (a(i+1),i=cc,cc+j-2,2)
c         write(6,500) "err", (sqrt(abs(covar(i+1,i+1))),i=cc,cc+j-2,2)
 42   continue

c      endif
 500  format(A3,1X,20(F11.6))
 501  format(A3,1X,F11.5,F11.5)
C      CALL PGCLOS()
C      rmavg=rmavg+avg

      GOTO 999

C      close(10)
C      goto 999

 901  write(6,*) 'could not open ', filename
      goto 999
 902  write(6,*) 'error in file called' , filename
      goto 999

 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca, 
     *     chisq,alamda,w,nstar) 
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,nca,ndata,ia(ma),MMAX 
      REAL alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca), 
     *     sig(ndata),x(ndata),y(ndata) 
      PARAMETER (MMAX=1000) 
      INTEGER j,k,l,mfit 
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX) 
      SAVE ochisq,atry,beta,da,mfit 
      if(alamda.lt.0.)then 
         mfit=0 
         do 11 j=1,ma 
            if (ia(j).ne.0) mfit=mfit+1 
 11      enddo  
         alamda=0.001 
         call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,w,
     .        nstar) 
         ochisq=chisq 
         do 12 j=1,ma 
            atry(j)=a(j) 
 12      enddo 
      endif 
      do 14 j=1,mfit 
         do 13 k=1,mfit
            covar(j,k)=alpha(j,k) 
 13      enddo 
         covar(j,j)=alpha(j,j)*(1.+alamda) 
         da(j)=beta(j) 
 14   enddo 
      call gaussj(covar,mfit,nca,da,1,1)  

      if(alamda.eq.0.)then 
         call covsrt(covar,nca,ma,ia,mfit) 
         call covsrt(alpha,nca,ma,ia,mfit) 
         return 
      endif 
      j=0 
      do 15 l=1,ma 
         if(ia(l).ne.0) then 
            j=j+1 
            atry(l)=a(l)+da(j) 
         endif 
 15   enddo 
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,w,nstar) 
      if(chisq.lt.ochisq)then 
         alamda=0.1*alamda 
         ochisq=chisq 
         do 17 j=1,mfit 
            do 16 k=1,mfit 
               alpha(j,k)=covar(j,k) 
 16         enddo 
            beta(j)=da(j) 
 17      enddo 
         do 18 l=1,ma 
            a(l)=atry(l) 
 18      enddo 
      else 
         alamda=10.*alamda 
         chisq=ochisq 
      endif 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp, 
     *     chisq,w,nstar) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,nalp,ndata,ia(ma),MMAX,nstar 
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata), 
     *     y(ndata),w
      PARAMETER (MMAX=1000) 
      INTEGER mfit,i,j,k,l,m 
      REAL dy,sig2i,wt,ymod,dyda(MMAX) 
      mfit=0 
      do 11 j=1,ma 
         if (ia(j).ne.0) mfit=mfit+1 
 11   enddo 
      do 13 j=1,mfit 
         do 12 k=1,j
            alpha(j,k)=0. 
 12      enddo 
         beta(j)=0. 
 13   enddo 
      chisq=0. 
      do 16 i=1,ndata 
         call funcs(x(i),a,ymod,dyda,ma,w,nstar) 
         sig2i=1./(sig(i)*sig(i)) 
         dy=y(i)-ymod 
         j=0 
         do 15 l=1,ma 
            if(ia(l).ne.0) then 
               j=j+1 
               wt=dyda(l)*sig2i 
               k=0 
               do 14 m=1,l 
                  if(ia(m).ne.0) then 
                     k=k+1 
                     alpha(j,k)=alpha(j,k)+wt*dyda(m) 
                  endif 
 14            enddo 
               beta(j)=beta(j)+dy*wt 
            endif 
 15      enddo 
         chisq=chisq+dy*dy*sig2i  
 16   enddo
      do 18 j=2,mfit 
         do 17 k=1,j-1 
            alpha(k,j)=alpha(j,k) 
 17      enddo
 18   enddo
      return 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE funcs(x,a,y,dyda,na,w,nstar) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER na,nstarmax
      parameter(nstarmax=100)
      REAL x,y,a(na),dyda(na) 
      real arg, cosarg,w(nstarmax),sinarg
      integer i,j,nstar,k,cc


      y=0
      j=(na-nstar)/nstar
      do 20 k=1,nstar
         cc=(j+1)*(k-1)+2
         dyda(cc-1)=0.
         do 21 i=cc,cc+j-2,2
            arg=real((i-cc+2)/2)*a(cc-1)*x+a(i+1)
            cosarg=cos(arg)
            sinarg=sin(arg)
            y=y+cosarg*a(i)
            dyda(i)=cosarg
            dyda(i+1)=-a(i)*sinarg
            dyda(cc-1)=dyda(cc-1)-a(i)*real((i-cc+2)/2)*x*sinarg
 21      continue
 20   continue

c      j=(na-nstar)/nstar
c      do 42 k=1,nstar
c         cc=(j+1)*(k-1)+2
c         write(6,500) "Amp", (a(i)  ,i=cc,cc+j-2,2)
c         write(6,500) "Psi", (a(i+1),i=cc,cc+j-2,2)
c 42   continue

 500  format(A3,1X,20(F9.6))

C     OLD STUFF FOR SINGLE FREQUENCY FITTING
C      y=0.
CC      y=a(1)
CC      dyda(1)=0.
C      do 11 i=1, NA-1, 2
C         arg=real((I+1)/2)*w*x+a(i+1)
C         cosarg=cos(arg)
C         y=y+cosarg*a(i)
C         dyda(i)=a(i)*cosarg
C         dyda(i+1)=-a(i)*sin(arg)
C 11   continue

      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER m,mp,n,np,NMAX 
      REAL a(np,np),b(np,mp) 
      PARAMETER (NMAX=1000) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), 
     *     ipiv(NMAX) 
      REAL big,dum,pivinv 
      do 11 j=1,n 
         ipiv(j)=0 
 11   enddo 
      do 22 i=1,n 
         big=0. 
         do 13 j=1,n 
            if(ipiv(j).ne.1)then 
               do 12 k=1,n 
                  if (ipiv(k).eq.0) then 
                     if (abs(a(j,k)).ge.big)then 
                        big=abs(a(j,k)) 
                        irow=j 
                        icol=k 
                     endif 
                  else if (ipiv(k).gt.1) then 
                     pause 'singular matrix in gaussj' 
                  endif 
 12            enddo 
            endif 
 13      enddo 
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then 
            do 14 l=1,n 
               dum=a(irow,l) 
               a(irow,l)=a(icol,l) 
               a(icol,l)=dum 
 14         enddo
            do 15 l=1,m 
               dum=b(irow,l) 
               b(irow,l)=b(icol,l) 
               b(icol,l)=dum 
 15         enddo 
         endif 
         indxr(i)=irow 
         indxc(i)=icol 
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj' 
         pivinv=1./a(icol,icol) 
         a(icol,icol)=1. 
         do 16 l=1,n 
            a(icol,l)=a(icol,l)*pivinv 
 16      enddo 
         do 17 l=1,m 
            b(icol,l)=b(icol,l)*pivinv 
 17      enddo 
         do 21 ll=1,n 
            if(ll.ne.icol)then 
               dum=a(ll,icol) 
               a(ll,icol)=0. 
               do 18 l=1,n 
                  a(ll,l)=a(ll,l)-a(icol,l)*dum 
 18            enddo
               do 19 l=1,m 
                  b(ll,l)=b(ll,l)-b(icol,l)*dum 
 19            enddo 
            endif 
 21      enddo
 22   enddo 
      do 24 l=n,1,-1 
         if(indxr(l).ne.indxc(l))then 
            do 23 k=1,n 
               dum=a(k,indxr(l)) 
               a(k,indxr(l))=a(k,indxc(l)) 
               a(k,indxc(l))=dum 
 23         enddo 
         endif 
 24   enddo 
      return 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,mfit,npc,ia(ma) 
      REAL covar(npc,npc) 
      INTEGER i,j,k 
      REAL swap
      do 12 i=mfit+1,ma 
         do 11 j=1,i 
            covar(i,j)=0. 
            covar(j,i)=0. 
 11      enddo 
 12   enddo 
      k=mfit 
      do 15 j=ma,1,-1 
         if(ia(j).ne.0)then 
            do 13 i=1,ma 
               swap=covar(i,k) 
               covar(i,k)=covar(i,j) 
               covar(i,j)=swap 
 13         enddo 
            do 14 i=1,ma 
               swap=covar(k,i) 
               covar(k,i)=covar(j,i) 
               covar(j,i)=swap 
 14         enddo 
            k=k-1 
         endif 
 15   enddo 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
         s=s+data(j)
 11   enddo
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         adev=adev+abs(s)
         p=s*s
         var=var+p
         p=p*s
         skew=skew+p
         p=p*s
         curt=curt+p
 12   enddo
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
         skew=skew/(n*sdev**3)
         curt=curt/(n*var**2)-3.
      else
         pause 'no skew or kurtosis when zero variance in moment'
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