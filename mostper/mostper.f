C23456789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7
      program mostper
      implicit none
C     MOST Fourier analysis software
C     (c) Jason Rowe 2007      
      integer nmax,npt,id,fatal,nb,nc,steps,nfitm,nfit,nfmax,idstart,m,
     .     bins,nfitl,nstar,nstarmax,plot,niter,i,ma,filet,iargc,fixfreq
      parameter(nmax=650000,nfitm=2000,nfmax=2000,nstarmax=300)
      real time(nmax),mag(nmax),merr(nmax),xcoo,ycoo,per1,per2,bper,
     .     a(nfmax),avg,fint,integrate,btheta,res(nmax),
     .     rmavg,bper1,bper2,btheta1,btheta2,tbin,w(nstarmax),
     .     pers(nstarmax),sig,amag(nmax),sky(nmax),xc(nmax),yc(nmax),
     .     fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),peak(nmax),nyquest,
     .     nyq,ofac,maxx,minx,aerr(nfmax),ftotal(nmax),aflux(nmax),
     .     nap(nmax),itime(nmax),mintime,skystd(nmax),sn,sns(nstarmax),
     .     snlimit
      double precision dtime(nmax),MOSTtime
      character*80 filename,ename,ans,cline
      logical multi
      common /data/ time,mag,merr,npt,nfit
      
      fixfreq=0 !0=fit frequencies, 1=fixed frequencies
      MOSTtime=2451545.00000000
      tbin=0.0
      snlimit=3.6
      per1=-1.0
      steps=0
      nfit=1  !for harmonic fits
      if(iargc().lt.2) goto 903
      call getarg(1,filename)
      call getarg(2,cline)
      read(cline,*) per2
      if(iargc().ge.3) then
         call getarg(3,cline)
         read(cline,*) per1
      endif
      if(iargc().ge.4) then
         call getarg(4,cline)
         read(cline,*) snlimit
      endif
      if(iargc().ge.5) then
         call getarg(5,cline)
         read(cline,*) tbin
      endif
      if(iargc().ge.6) then
         call getarg(6,cline)
         read(cline,*) steps
      endif
      if(iargc().ge.7) then
         call getarg(7,cline)
         read(cline,*) nfit
      endif
      if(iargc().ge.8) then
         call getarg(8,cline)
         read(cline,*) fixfreq
      endif
      if(nfit.lt.1) goto 903
      if((fixfreq.lt.0).or.(fixfreq.gt.2)) goto 904
      multi=.false.
      rmavg=0.0

c      write(6,*) "Enter name of filename"
c      read(5,500) filename
 500  format(a80)
      open(unit=10,file=filename,status='old',err=901)

      open(unit=14,file="freqs.dat")

C     filet=0 for Direct, filet=1 for Fabry
      filet=0
      call readdata(fatal,npt,id,xcoo,ycoo,time,dtime,mag,merr,sky,xc,
     .   yc,fx,fy,fxy,tboard,peak,ftotal,aflux,nap,itime,skystd,filet,
     .   mintime,MOSTtime)
     
C     These lines convert to mmag
      do 11 i=1,npt
      	mag(i)=mag(i)*1000.0
        merr(i)=merr(i)*1000.0
 11   continue     
     
      write(6,501) "Zero Time: ",mintime
      write(6,*) "Points Read: ",npt
 501  format(A11,F13.8)
      id=1

      minx= 99.9e30
      maxx=-99.9e30
      do 5 i=1,npt
         maxx=max(time(i),maxx)
         minx=min(time(i),minx)
 5    continue

      call pgopen('?')
c      call pgask(.false.)
      call pgsch(2.9)
      call pgsubp(1,4)
      call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels

      call pgpage()


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Commands Start Here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     plotdata
C     pdm
C     plotph
C     cosinefit
C     grelor (not yet functional)
C     windowfn
C     avgrm
C     shaperm

      nyq=nyquest(npt,time)

      nb=5
      nc=2
      if(per1.lt.0.0) per1=1.0/(maxx-minx)
c      per1=3.0
c      per2=20.0
c      per1=200.0
c      per2=250.0
      ofac=4.0
      write(6,*) "nyq",nyq
      if(steps.le.0) steps=int(ofac*(per2-per1)*npt/nyq)
c      steps=12000
      if(steps.gt.1e6)then
         write(6,*) "reducing steps to 1e6"
         steps=1e6
      endif
c      steps=45000
c      steps=4000
      write(6,*) "freq1,freq2,steps",per1,per2,steps
      bins=150
c      tbin=3.0*1440.0
      nfitl=3
      niter=nstarmax-1
      sig=3.0

      do 10 i=1,0
      call sigclip(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,fxy,
     .     tboard,peak,ftotal,aflux,nap,itime,rmavg)
 10   continue

c      call bind(npt,time,mag,merr,bins)
      call avgrm(npt,mag,rmavg)
      if(tbin.gt.0.0) call bindt(npt,time,mag,merr,tbin)
c      call detrend(npt,time,mag,merr,nfitl)
c      tbin=1.0*500.0
      tbin=0.25*24.0*60.0
c      call spdetrend(npt,time,mag,merr,tbin,sky)
c      call sigclip(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,fxy,
c     .     tboard,peak,ftotal,aflux,nap,itime,rmavg)
c      goto 21
c      call plotdata(npt,time,mag,merr,-1.0,-1.0,1,1,MOSTtime)
c      write(6,*) id,xcoo,ycoo,npt
      call avgrm(npt,mag,rmavg)

      write(0,*) "rmavg: ",rmavg

c      do 100 i=1,npt
c         itime(i)=tbin*60.0
c 100  continue

      nstar=0
      
c      goto 23

c      call plotph(npt,time,mag,3.52474859,1,2)
      call plotdata(npt,time,mag,merr,-1.0,-1.0,1,1,MOSTtime)
c      goto 23
      call jmfourw(npt,time,mag,merr,per1,per2,steps,bper,btheta,sn,1,2,
     ,   snlimit)
c      goto 23
c      call pdm(npt,time,mag,nb,nc,per1,per2,steps,bper,btheta,1,1)
c      bper=1.0/bper
c      call plotph(npt,time,mag,bper,1,2)
c      write(6,*) "bper:",bper!,sn
c      goto 23
c      bper=1.0/0.521183431
c      bper=1.0/0.0704256
c      call pdm(npt,time,mag,nb,nc,1.8,2.0,2000,bper,btheta,1,3)
c      bper=1.0/12.931
      write(6,*) "bper:",bper
      plot=0
      bper=1.0/bper
c      call plotph(npt,time,mag,bper,1,3)
c      goto 900
c      bper=3.52474859
      nstar=nstar+1
      pers(nstar)=bper
      sns(nstar)=sn
      call cosinefit(pers,res,id,xcoo,ycoo,avg,a,btheta,rmavg,plot,
     .     ma,aerr,1,3,w,nstar,fixfreq)
      call windowfn(npt,nfit,time,mag,merr,bper,per1,per2,avg,a,nb,nc,
     .   steps,1,3,nstar,pers)
      plot=0
      call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,
     .     nstar,MOSTtime)

c       call plotph(npt,time,mag,bper,1,4)
c      per2=100.0
c      steps=int(ofac*(per2-per1)*npt/nyq)
      call jmfourw(npt,time,res,merr,per1,per2,steps,bper,btheta,sn,1,4,
     .  snlimit)

c      call plotph(npt,time,mag,bper,1,3)


c      call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
c      call windowfn(npt,nfit,time,mag,bper,per1,per2,avg,a,nb,nc,steps,
c     .     1,4)
c      call jmfourw(npt,time,res,merr,per1,per2,steps,bper,btheta,1,4,
c     .  snlimit)
c      call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,1,4)

      plot=1
      write(6,*) "Next peak?          "
      read(5,*) ans
      if(ans.eq."n")goto 24

      i=0

      niter=40
  20  i=i+1
      if((sn.lt.snlimit).or.(i.gt.niter)) then
         ans="n"
      else
         ans="y"
      endif
c      if(i.gt.niter) then
c         ans="n"
c      else
c         ans="y"
c      endif

      if(ans.eq."y")then
         write(6,*) "Iteration #: ",i
         write(6,*) "bfreq:",bper,sn
cp         call pgpanl(1,4)
cp         call pgpage()
         bper=1.0/bper
         nstar=nstar+1
         pers(nstar)=bper
         sns(nstar)=sn
cp         call jmfour(npt,time,mag,per1,per2,steps,bper,btheta,1,1)
         call pgpanl(1,2)
         call pgeras()
c         call plotph(npt,time,res,bper,1,2)
         plot=1
         call cosinefit(pers,res,id,xcoo,ycoo,avg,a,btheta,rmavg,plot,
     .        ma,aerr,1,2,w,nstar,fixfreq)
         write(6,503) i+1,"bfreq:",1.0/bper,sn
 503     format(i3,1x,A6,1X,F8.5,1X,F6.3)
c         call plotdata(npt,time,mag,merr,-1.0,-1.0,1,2,MOSTtime)
c         call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
         call pgpanl(1,3)
         call pgeras()
         call shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,1,3,plot,w,
     .        nstar,MOSTtime)
         call pgpanl(1,4)
c         call pgeras()
         call jmfourw(npt,time,res,merr,per1,per2,steps,bper,btheta,sn,
     .      1,4,snlimit)
c         call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,1,4)
         goto 20
      endif
      
C     OLDIES BUT GOODIES
c      call pdm(npt,time,mag,nb,nc,per1,per2,steps,bper,btheta,1,2)
C      m=12
c      call grelor(npt,time,mag,merr,m,per1,per2,steps,1,3)
c      call windowfn(npt,nfit,time,mag,bper,per1,per2,avg,a,nb,nc,steps,
c     .     1,4)


 24   continue
 
      if (fixfreq.eq.2) call updatecosinefit(a,rmavg,ma,nstar,aerr) 
      
c      call pgpanl(1,1)
c      call pgeras()
c      call plotdata(npt,time,mag,merr,-1.0,-1.0,1,1,MOSTtime)
c      call plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
c      call bindt(npt,time,res,merr,tbin)
      
 21   write(14,502) mintime
 502  format(F13.8)
      call freqout(ma,nstar,a,aerr,sns)
      close(14)
      write(6,501) "Zero Time: ",mintime
 23   ename="test.dat"   

      if(tbin.gt.0.0)then
         rmavg=0.0
         do 22 i=1,npt
            dtime(i)=time(i)+mintime
 22      continue
      endif
   
      write(6,*) "rmavg:",rmavg
C     These lines convert to mmag
      do 12 i=1,npt
        mag(i)=mag(i)/1000.0
        merr(i)=merr(i)/1000.0
 12   continue    
      call exportdata(npt,dtime,mag,merr,ename,sky,xc,yc,
     .      fx,fy,fxy,tboard,peak,ftotal,aflux,nap,itime,skystd,rmavg,
     .      filet,mintime)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      call pgpanl(1,4)
c      read(5,*)
c      call pgpage()
c      goto 10


 900  call pgclos()
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 902  close(10)
      call pgclos()
      write(6,*) "End of datafile reached..."
      goto 999
 903  write(6,*) "Usage: mostper <filename> <per2> [per1] [sn] [tbin] [s
     .teps] [nfit] [fixfreq]"
      write(6,*) "<filename>: MOST photometry file"
      write(6,*) "<per2>: largest frequency to scan"
      write(6,*) "[per1]: lowest frequency to scan (optional)"
      write(6,*) "        set per1 negative for default"
      write(6,*) "[sn]: S/N for statistically significant frequencies (o
     .ptional)"
      write(6,*) "        default = 3.6"
      write(6,*) "[tbin]: Binning time scale (minutes) (optional)"
      write(6,*) "        set to < 0 to disable binning"
      write(6,*) "[steps]: Steps for DFT scan"
      write(6,*) "         default is 4 times oversampling"
      write(6,*) "[nfit]: number of harmonics to fit (optional)"
      write(6,*) "        default is 1"
      write(6,*) "[fixfreq]: 0 (default): fit frequencies"
      write(6,*) "           1 : do not fit frequencies"
      write(6,*) "           2 : do not fix frequencies, except for fina
     .l iteration"
      goto 999
 904  write(6,*) "Bad Usage:"
      write(6,*) "[fixfreq] must be 0, 1 or 2."
      goto 999
 999  end


Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(fatal,npt,id,xcoo,ycoo,time,dtime,mag,merr,
     .   sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,
     .   skystd,filet,mintime,MOSTtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     id - id of star (returned)
C     xcoo,ycoo - co-ordinates of star (returned)
C     time,mag,merr - the data and associated error.
      implicit none
      integer nmax,nce
      parameter(nmax=650000,nce=11)
      integer fatal,npt,id,i,filet,rn,j
      real xcoo,ycoo,time(nmax),mag(nmax),merr(nmax),c,amag(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),skystd(nmax),mintime,ntest
      double precision dtime(nmax),cevents(nce),twid,MOSTtime,dmintime
      logical flag
      data cevents /2045.67,2047.78,2048.84,2050.86,2056.99,2059.17,
     .   2061.35,2063.67,2075.40,2077.20,2081.87/
   
      fatal=0

C      read(10,*,end=30) id,xcoo,ycoo

      i=1
      twid=0.25

      dmintime=99.9d30
      if(filet.eq.0) then
   9     continue
  10     read(10,*,err=9,end=20) dtime(i),mag(i),merr(i),sky(i),
     .      xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .      ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
c             merr(i)=0.0001
c            mag(i)=sky(i)
            time(i)=dtime(i)
            dmintime=min(dtime(i),dmintime)
            
            flag=.false.
c            do 12 j=1,nce
c               if((dtime(i).gt.cevents(j)-twid).and.
c     .           (dtime(i).lt.cevents(j)+twid)) flag=.true.
c 12         continue            
c            if(flag)goto 10            
            
c            if(mfield(i).lt.20000.0) goto 10
c            if(itime(i).lt.16.0) goto 10
c            merr(i)=0.001
c            if(mag(i).gt.0.05) goto 10
c            if(mag(i).lt.-0.05) goto 10
c            if((sky(i).lt.0.0).or.(sky(i).gt.2000.0)) goto 10
            merr(i)=abs(merr(i))
            if(merr(i).le.0.0) goto 10
c            if(sky(i).gt.1500.0) goto 10
c         if((time(i).gt.1690.4).and.(time(i).lt.1690.6)) goto 10
c         if((time(i).gt.1693.8).and.(time(i).lt.1694.2)) goto 10
c         if((time(i).gt.1697.4).and.(time(i).lt.1697.7)) goto 10
c         if((time(i).gt.1700.9).and.(time(i).lt.1701.2)) goto 10
         rn=ntest(mag(i))
         if(rn.eq.1) goto 10
         rn=ntest(sky(i))
         if(rn.eq.1) goto 10
            i=i+1
         goto 10
 20      continue
      elseif(filet.eq.1) then
 29      continue
 30      read(10,*,err=29,end=40) dtime(i),mag(i),merr(i),sky(i)
C     .      xc(i),yc(i)!,fx(i),fy(i)
           time(i)=dtime(i)
           mintime=min(time(i),mintime)
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
      
      mintime=real(dmintime)
      MOSTtime=MOSTtime+dmintime !calculate real JD offset
      do 50 i=1,npt
         time(i)=time(i)-mintime
 50   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine freqout(ma,nstar,a,aerr,sns)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ma,nstar,nfmax,i,j,k,cc,b,nstarmax
      parameter(nfmax=2000,nstarmax=300)
      real a(nfmax),twopi,pi,aerr(nfmax),sns(nstarmax)

      pi=3.141592654
      twopi=2.0*pi

      j=(ma-nstar)/nstar
      do 38 k=1,nstar
         cc=(j+1)*(k-1)+2
         do 39 i=cc,cc+j-2,2
            if(a(i).lt.0.0) then
               a(i)=abs(a(i))
               a(i+1)=a(i+1)+pi
            endif
            b = int(a(i+1)/twopi)
            a(i+1)=a(i+1)-b*twopi
            if (a(i+1).lt.0) a(i+1)=twopi+a(i+1)
 39      continue
 38   continue

      write(14,502) ma,nstar
      
      j=(ma-nstar)/nstar
      do 10 k=1,nstar
         cc=(j+1)*(k-1)+2
         write(6,500) a(cc-1)/twopi,(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
         write(14,500) a(cc-1)/twopi,(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
         write(6,500) aerr(cc-1),(aerr(i),i=cc,cc+j-2,2),
     .        (aerr(i+1),i=cc,cc+j-2,2)
 10   continue

 500  format(20(E14.7,1X))
c 500  format(20(F13.6,1X))
 501  format(A3,1X,F11.4)
 502  format(I4,1X,I4)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigclip(npt,time,mag,merr,sig,sky,xc,yc,fx,fy,
     .     fxy,tboard,peak,ftotal,aflux,nap,itime,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=650000,nitmax=50)
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
      subroutine plotfit(npt,nfit,time,mag,pers,avg,a,w,nstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit,nstarmax,nstar,i,j,k,i2,cc,nmax,na,steps
      parameter(nstarmax=300,nmax=650000,steps=2000)
      real time(npt),mag(npt),pers(nstarmax),avg,a(nfit),w(nstarmax),pi,
     .     temp,px(nmax),py(nmax),arg,tt,maxt,mint

      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar

      maxt=-99.9e30
      mint= 99.9e30
      do 7 i=1,npt
         maxt=max(maxt,time(i))
         mint=min(mint,time(i))
 7    continue

      j=(na-nstar)/nstar
      do 10 i2=1,steps
         tt=real(i2)*(maxt-mint)/real(steps)+mint
         temp=avg
c         write(6,*) avg
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
c               write(6,*) cc-1,(i-cc+2)/2,a(cc-1),a(i),a(i+1)
               arg=real((i-cc+2)/2)*a(cc-1)*tt+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         px(i2)=tt
         py(i2)=temp
 10   continue
      call pgsci(2)
      call pgline(steps,px,py)
      call pgsci(1)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .      fx,fy,fxy,tboard,peak,ftotal,aflux,nap,itime,skystd,rmavg,
     .      filet,mintime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,filet
      real time(npt),mag(npt),merr(npt),sky(npt),xc(npt),
     .   yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),peak(npt),rmavg,
     .   ftotal(npt),aflux(npt),nap(npt),itime(npt),mintime,
     .   skystd(npt)
      double precision dtime(npt)
      character*80 filename

      open(unit=11,file=filename)

      do 10 i=1,npt
         if(filet.eq.1) then
            xc(i)=0.
            yc(i)=0.
            fx(i)=0.
            fy(i)=0.
            fxy(i)=0.
            tboard(i)=0.
            peak(i)=0.
            aflux(i)=0.
            nap(i)=0.
            itime(i)=0.
         endif
         write(11,510) dtime(i),mag(i)+rmavg,merr(i),sky(i),
     .        xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),peak(i),
     .        ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
 10   continue

 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)

      close(11)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine spdetrend(npt,time,mag,merr,tbin,sky)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i,bins,bnpt,sbins
      parameter(nmax=650000)
      real time(npt),mag(npt),merr(npt),bdata(nmax),bderr(nmax),yp1,
     .     yp2,y2(nmax),btime(nmax),x,y,tbin,sky(npt)

      bins=sbins

      bnpt=0
      do 10 i=1,npt
c         if(sky(i).lt.2000.0) then
            bnpt=bnpt+1
            btime(bnpt)=time(i)
            bdata(bnpt)=mag(i)
            bderr(bnpt)=merr(i)
c         endif
 10   continue

      call bindt(bnpt,btime,bdata,bderr,tbin)

      yp1=1.0e30
      yp2=1.0e30
      call spline(btime,bdata,bnpt,yp1,yp2,y2)
      
      do 20 i=1,npt
         x=time(i)
c         write(6,*) x,y
         call splint(btime,bdata,y2,bnpt,x,y)
         mag(i)=mag(i)-y
 20   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine detrend(npt,time,mag,merr,nfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit
      integer ia(nfit),i,ma,j
      real time(npt),mag(npt),merr(npt),ans(nfit),covar(nfit,nfit),
     .     chisq,T

      do 36 i=1,nfit
         ia(i)=1
 36   continue

      ma=nfit
      call lfit(time,mag,merr,npt,ans,ia,ma,covar,nfit,chisq)
      write(6,*) (ans(i),i=1,nfit)

      do 10 i=1,npt
         T=0.
         DO 33 J = 1, NFIT
            T = ANS(J)*((time(i))**REAL(J-1)) +  T
 33      CONTINUE
         mag(i)=mag(i)-T
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fourcal(npt,time,mag,per1,per2,steps,bper,btheta,
     .     panx,pany,w,nstar,rmavg,fixfreq)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,i,nstarmax,n,nfit,nstar,id,
     .     plot,ma,nfitm,nfmax,j,k,cc,l,jj,stepmax,fixfreq
      parameter(nmax=650000,nstarmax=300,nfitm=2000,nfmax=2000,
     .     stepmax=1000)
      real time(npt),mag(npt),per1,per2,bper,btheta,ptest,w(nstarmax),
     .     x(nmax),y(nmax),sig(nmax),pers(nstarmax),
     .     res(nmax),avg,a(nfmax),rmavg,aerr(nfmax),twopi,pi,px(3),
     .     py(3),amp,amperr,phi,phierr,nx(stepmax),ny1(stepmax),
     .     ny2(stepmax),perc,percold,sn,xcoo,ycoo
      character*80 tmpc
      common /data/ x,y,sig,n,nfit

      pi=3.141592654
      twopi=2.0*pi
      ptest=per1

      id=1
      xcoo=1.0
      ycoo=1.0
      nfit=1
      avg=0

      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

      plot=0
      btheta=-99.9e10

c      open(unit=12,file="t1.dat")

      do 10 i=1,steps

         perc=100.0*real(i)/real(steps)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
 500        format(A13,F6.1)
            percold=perc
            call ovrwrt(tmpc,2)
         endif
         
         jj=1+real(stepmax)/real(steps)*(i-1)

         nstar=1
         pers(1)=1.0/ptest
         call cosinefit(pers,res,id,xcoo,ycoo,avg,a,btheta,rmavg,
     .        plot,ma,aerr,panx,pany,w,nstar,fixfreq)

         call freqout(ma,nstar,a,aerr,sn)
         j=(ma-nstar)/nstar
         do 42 k=1,nstar
            cc=(j+1)*(k-1)+2
            amp=a(cc)
            amperr=aerr(cc)
            phi=a(cc+1)
            phierr=aerr(cc+1)
            write(6,*) jj,a(cc-1)/twopi,(a(l),aerr(l),a(l+1),aerr(l+1),
     .           l=cc,cc+j-2,2)
 42      continue
         ptest=ptest+(per2-per1)/real(steps)
         
c         write(12,*) a(cc-1)/twopi,amp

         nx(jj)=ptest
         if(amp.gt.ny1(jj)) ny1(jj)=amp
         if(amp.lt.ny2(jj)) ny2(jj)=amp
         if(amp.gt.btheta) then
            btheta=amp
            bper=ptest
         endif
 10   continue
      
c      close(12)

      call pgpanl(panx,pany)
      call pgwindow(per1,per2,0.0,btheta)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgbbuf()
      do 30 i=2,stepmax
         px(1)=nx(i-1)
         py(1)=ny2(i-1)
         px(2)=nx(i)
         px(3)=nx(i)
         py(2)=ny1(i)
         py(3)=ny2(i)
         call pgline(3,px,py)
 30   continue
      call pgebuf()

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jmfourw(npt,time,mag,merr,per1,per2,steps,bper,btheta,
     .     signoise,panx,pany,snlimit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,stepmax,k,j,nf,snw,n1,n2,ndb
      parameter(stepmax=2000,snw=2000,ndb=1000000)
      real time(npt),mag(npt),per1,per2,bper,btheta,percold,perc,by(4),
     .     nx(stepmax),ny1(stepmax),ny2(stepmax),px(3),py(3),bx(4),tbin,
     .     orbper,dum,weight,merr(npt),amps(ndb),cf(2*snw+2),stdev,std,
     .     cf2(2*snw+2),mean,signoise,snplot(ndb),freqs(ndb),snlimit,
     .     dumr
      parameter(nmax=650000)
      character*80 tmpc

      INTEGER I,NUMPTS,NMDENU
      DOUBLE PRECISION JD (nmax),B(nmax),LONU,DELTNU,TWOPI
      DOUBLE PRECISION SIDEL(nmax),CODEL(nmax),SI(nmax),CO(nmax)
      DOUBLE PRECISION SISUM,COSUM,SN,NU,AMP,PHI
      CHARACTER *20 INFILE,SPFILE
      DOUBLE PRECISION JDTMP,LOJD,HIJD,BTMP,JDZRO

c      snlimit=3.6

      if(steps.lt.stepmax) then
         write(6,*) "Increasing number of steps to 2000"
         steps=stepmax
      endif


      tbin=0.0

      btheta=-99.9e10
      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

      open(unit=21,file="ft.dat")

      LOJD=0.
      HIJD=100.
      JDZRO=0.

      LONU=dble(per1)
      DELTNU=dble(per2-per1)/dble(steps)
      NMDENU=steps

      do 3 i=1,npt
         jd(i)=dble(time(i))
         b(i)=dble(mag(i))
 3    continue
      NUMPTS = npt

      TWOPI = 6.2831853
      NU = LONU
      DO 61 I=1,NUMPTS
        SI(I) = SIN(TWOPI*LONU*JD(I))
        CO(I) = COS(TWOPI*LONU*JD(I))
        SIDEL(I) = SIN(TWOPI*DELTNU*JD(I))
        CODEL(I) = COS(TWOPI*DELTNU*JD(I))
   61 CONTINUE

      percold=0.0      
      DO 51 J=1,NMDENU
         perc=100.0*real(j)/real(steps)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
            percold=perc
            call ovrwrt(tmpc,2)
         endif
 500     format(A13,F6.1)
         SISUM = 0.
         COSUM = 0.
         weight=0.0
         DO 41 I=1,NUMPTS
            SISUM = SISUM + B(I)*SI(I)/merr(i)
            COSUM = COSUM + B(I)*CO(I)/merr(i)
            weight=weight+1.0/merr(i)
 41      CONTINUE
C     WRITE (50,*) NU,SISUM,COSUM
         AMP = (2./weight)*SQRT(SISUM**2. + COSUM**2.)
         PHI = ATAN(-SISUM/COSUM)
C     Jaymie is stupid
c         PHI= ATAN(COSUM/SISUM)
         IF (COSUM.LT.0) PHI = PHI + TWOPI/2.
         IF ((SISUM.GT.0).AND.(COSUM.GE.0)) PHI = PHI + TWOPI
         amps(j)=real(amp)
         freqs(j)=real(nu)
         k=1+real(stepmax)/real(steps)*(j-1)
         nx(k)=real(nu)
         if(real(amp).gt.ny1(k)) ny1(k)=real(amp)
         if(real(amp).lt.ny2(k)) ny2(k)=real(amp)
         if(tbin.gt.0.) then
         if((real(amp).gt.btheta).and.(real(nu).gt.1440.0/tbin)) then
            btheta=real(amp)
            bper=real(nu)
            nf=j
         endif
         else
         if(real(amp).gt.btheta) then
            btheta=real(amp)
            bper=real(nu)
            nf=j
         endif
         endif 
c     call pgpt1(NU,AMP,1)
      WRITE (21,*) NU,AMP,PHI
c     
         DO 31 I=1,NUMPTS
            SN = SI(I)
            SI(I) = SN*CODEL(I) + CO(I)*SIDEL(I)
            CO(I) = CO(I)*CODEL(I) - SN*SIDEL(I)
 31      CONTINUE
         NU = NU + DELTNU
 51   CONTINUE
 
cC     do s/n calc.
c      n1=nf-snw-1
c      n2=nf+snw
c      if(n1.le.0)n1=1
c      if(n2.gt.NMDENU)n2=NMDENU
c      j=0
c      do 200 i=n1,n2
c         j=j+1
c         cf(j)=amps(i)
c 200  continue
c      std=stdev(j,cf,mean)
c      k=0
c      do 201 i=1,j
c         if(abs(cf(i)-mean).lt.3.0*std) then
c            k=k+1
c            cf2(k)=cf(i)
c         endif
c 201  continue
c      std=stdev(k,cf2,mean)
c      signoise=btheta/mean
      
c      call pgpage()
      call pgpanl(panx,pany)
      call pgeras()
      call pgwindow(per1,per2,0.0,btheta)
cc      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pgbox('BTS',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Frequency","AMP"," ")
      call pgptxt((per1+per2)/2.0,-0.25*(btheta),0.0,0.5,
     .	"Frequency (c/d)")
      call pgptxt((per1+per2)/2.0,btheta+0.23*(btheta),0.0,0.5,
     .	"Frequency (mHz)")
      call pgptxt(per1-0.05*(per2-per1),(btheta)/2,90.0,0.5,
     .	"Amplitude (\(0638)mag)")
      call pgbbuf()
      
      close(21)
      do 30 i=2,stepmax
         px(1)=nx(i-1)
         py(1)=ny2(i-1)
         px(2)=nx(i)
         px(3)=nx(i)
         py(2)=ny1(i)
         py(3)=ny2(i)
         call pgline(3,px,py)
 30   continue
      call pgebuf()


      call pgsls(2)
      call pgsci(2)
      orbper=14.199363
      i=per1/orbper+1
      j=per2/orbper
      
      do 100 k=i,j
         px(1)=k*orbper
         px(2)=px(1)
         py(1)=0.0
         py(2)=btheta
         call pgline(2,px,py)
 100  continue

      call pgsls(1)
      call pgsci(1)
      
      signoise=0.0
      btheta=0.0
C     do s/n plot
      do 302 nf=1,steps
         n1=nf-snw-1
         n2=nf+snw
         if(n1.le.0)n1=1
         if(n2.gt.steps)n2=steps
         j=0
         do 300 i=n1,n2
            j=j+1
            cf(j)=amps(i)
 300     continue
         std=stdev(j,cf,mean)
         k=0
         do 301 i=1,j
            if(abs(cf(i)-mean).lt.3.0*std) then
               k=k+1
               cf2(k)=cf(i)
            endif
 301     continue
         std=stdev(k,cf2,mean)
         snplot(nf)=snlimit*mean
         if(amps(nf)/mean.gt.snlimit)then  !is S/N good?
            if(amps(nf).gt.btheta)then  !pick largest amplitude
               signoise=amps(nf)/mean
               bper=freqs(nf)
               btheta=amps(nf)
            endif
         endif
c         write(6,*) n1,n2,nf,freqs(nf),mean
c         read(5,*)
 302  continue
      call pgsci(3)
c      call pgsls(2)
C     draw S/N line
      call pgline(steps,freqs,snplot)
      call pgsci(1)
c      call pgsls(1)

      call pgwindow(per1,per2,0.0,1000.0*btheta)
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)

      dumr=1000.0/86400.0
      call pgwindow(dumr*per1,dumr*per2,0.0,1000.0*btheta)
      call pgbox('CMTS1',0.0,0,'',0.0,0)
c      call pgbox('CTS',0.0,0,'',0.0,0)

c      px(1)=360.0/0.997268
c      px(2)=px(1)
c      py(1)=0.0
c      py(2)=btheta
c      call pgline(2,px,py)
c      call pgsls(1)
c      if(tbin.gt.0.0) then
c         bx(1)=per1
c         by(1)=0.0
c         bx(2)=(60.0*24.0)/tbin
c         by(2)=0.0
c         bx(3)=bx(2)
c         by(3)=btheta
c         bx(4)=per1
c         by(4)=btheta
c         call pgsfs(3)
c         call pgpoly(4,bx,by)
c         call pgsfs(1)
c      endif
c      call pgsci(1)
      
      return
      end

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine dft(npt,time,mag,per1,per2,steps,bper,btheta,panx,pany)
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      integer npt,stepmax,nmax,steps,panx,pany
c      parameter(stepmax=5000,nmax=650000)
c      real time(npt),mag(npt),per1,per2,bper,btheta      
c
c
c
c DO R = 1, FSZ
c   NUS(R) = FSTART + DBLE(R-1)*FSTEPB
c  END DO
c   TWZERO = TW_PI * FSTART
c   !CALCULATE SPECTRUM
c   SISUM=0.0D0
c   COSUM=0.0D0
c DO R = 1, TSZ
c   TO = TWZERO*TIME(R)
c   TO = DMOD(TO,TW_PI)
c   TD = TWSTEP*TIME(R)
c   TD = DMOD(TD,TW_PI)
c   SIDEL = DSIN(TD)
c   CODEL = DCOS(TD)
c   SOLD  = COUNTS(R)*DSIN(TO)
c   COLD  = COUNTS(R)*DCOS(TO)
c   SNEW  = SOLD
c   CNEW  = COLD
c   DO J = 1, FSZ
c    SISUM(J) = SISUM(J) + SNEW
c    COSUM(J) = COSUM(J) + CNEW
c    SNEW = SOLD * CODEL + COLD * SIDEL
c    CNEW = COLD * CODEL - SOLD * SIDEL
c    SOLD = SNEW
c    COLD = CNEW
c   ENDDO
c END DO
c
c  AMPS(:)   = (2.0D0/DBLE(TSZ))*DSQRT(SISUM*SISUM + COSUM*COSUM)
c  PHASES(:) = DATAN2 (COSUM,SISUM)
c  
c      return
c      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jmfour(npt,time,mag,per1,per2,steps,bper,btheta,
     .     panx,pany)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,steps,panx,pany,nmax,stepmax,k,j
      parameter(stepmax=5000)
      real time(npt),mag(npt),per1,per2,bper,btheta,percold,perc,by(4),
     .     nx(stepmax),ny1(stepmax),ny2(stepmax),px(3),py(3),bx(4),tbin,
     .     orbper,dum
      parameter(nmax=650000)
      character*80 tmpc

      INTEGER I,NUMPTS,NMDENU
      DOUBLE PRECISION JD (nmax),B(nmax),LONU,DELTNU,TWOPI
      DOUBLE PRECISION SIDEL(nmax),CODEL(nmax),SI(nmax),CO(nmax)
      DOUBLE PRECISION SISUM,COSUM,SN,NU,AMP,PHI
      CHARACTER *20 INFILE,SPFILE
      DOUBLE PRECISION JDTMP,LOJD,HIJD,BTMP,JDZRO

      if(steps.lt.stepmax) then
         write(6,*) "Increasing number of steps to 2000"
         steps=stepmax
      endif


      tbin=0.0

      btheta=-99.9e10
      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

      open(unit=21,file="ft.dat")

      LOJD=0.
      HIJD=100.
      JDZRO=0.

      LONU=dble(per1)
      DELTNU=dble(per2-per1)/dble(steps)
      NMDENU=dble(steps)

      do 3 i=1,npt
         jd(i)=dble(time(i))
         b(i)=dble(mag(i))
 3    continue
      NUMPTS = npt

      TWOPI = 6.2831853
      NU = LONU
      DO 61 I=1,NUMPTS
        SI(I) = SIN(TWOPI*LONU*JD(I))
        CO(I) = COS(TWOPI*LONU*JD(I))
        SIDEL(I) = SIN(TWOPI*DELTNU*JD(I))
        CODEL(I) = COS(TWOPI*DELTNU*JD(I))
   61 CONTINUE

      percold=0.0      
      DO 51 J=1,NMDENU
         perc=100.0*real(j)/real(steps)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
            percold=perc
            call ovrwrt(tmpc,2)
         endif
 500     format(A13,F6.1)
         SISUM = 0.
         COSUM = 0.
         DO 41 I=1,NUMPTS
            SISUM = SISUM + B(I)*SI(I)
            COSUM = COSUM + B(I)*CO(I)
 41      CONTINUE
c         WRITE (50,*) NU,SISUM,COSUM
         AMP = (2./NUMPTS)*SQRT(SISUM**2. + COSUM**2.)
         PHI = ATAN(-SISUM/COSUM)
C        Jaymie is stoopid
C         PHI = ATAN(COSUM/SISUM)
         IF (COSUM.LT.0) PHI = PHI + TWOPI/2.
         IF ((SISUM.GT.0).AND.(COSUM.GE.0)) PHI = PHI + TWOPI
c     
         k=1+real(stepmax)/real(steps)*(j-1)
         nx(k)=real(nu)
         if(real(amp).gt.ny1(k)) ny1(k)=real(amp)
         if(real(amp).lt.ny2(k)) ny2(k)=real(amp)
         if(tbin.gt.0.) then
         if((real(amp).gt.btheta).and.(real(nu).gt.1440.0/tbin)) then
            btheta=real(amp)
            bper=real(nu)
         endif
         else
         if(real(amp).gt.btheta) then
            btheta=real(amp)
            bper=real(nu)
         endif
         endif 
c     call pgpt1(NU,AMP,1)
      WRITE (21,*) NU,AMP,PHI
c     
         DO 31 I=1,NUMPTS
            SN = SI(I)
            SI(I) = SN*CODEL(I) + CO(I)*SIDEL(I)
            CO(I) = CO(I)*CODEL(I) - SN*SIDEL(I)
 31      CONTINUE
         NU = NU + DELTNU
 51   CONTINUE
      
c      call pgpage()
      call pgpanl(panx,pany)
      call pgwindow(per1,per2,0.0,btheta)
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Frequency","AMP"," ")
      call pgbbuf()
      
      close(21)
      do 30 i=2,stepmax
         px(1)=nx(i-1)
         py(1)=ny2(i-1)
         px(2)=nx(i)
         px(3)=nx(i)
         py(2)=ny1(i)
         py(3)=ny2(i)
         call pgline(3,px,py)
 30   continue
      call pgebuf()


      call pgsls(2)
      call pgsci(2)
      orbper=14.199363
      i=per1/orbper+1
      j=per2/orbper
      
      do 100 k=i,j
         px(1)=k*orbper
         px(2)=px(1)
         py(1)=0.0
         py(2)=btheta
         call pgline(2,px,py)
 100  continue

      call pgsls(1)
      call pgsci(1)

      dum=1000.0/86400.0
      call pgwindow(dum*per1,dum*per2,0.0,btheta)
      call pgbox('CMTS1',0.0,0,'',0.0,0)

c      px(1)=360.0/0.997268
c      px(2)=px(1)
c      py(1)=0.0
c      py(2)=btheta
c      call pgline(2,px,py)
c      call pgsls(1)
c      if(tbin.gt.0.0) then
c         bx(1)=per1
c         by(1)=0.0
c         bx(2)=(60.0*24.0)/tbin
c         by(2)=0.0
c         bx(3)=bx(2)
c         by(3)=btheta
c         bx(4)=per1
c         by(4)=btheta
c         call pgsfs(3)
c         call pgpoly(4,bx,by)
c         call pgsfs(1)
c      endif
c      call pgsci(1)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCS2(X,P,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Taken from NUMERICAL RECIPES: The Art of Scientific Computing by
C     Press, Flannery, Teukolsky, Vetterling (1986).
C
      DIMENSION P(NP)
      P(1)=1.
      DO 11 J=2,NP
         P(J)=P(J-1)*X
 11   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindt(npt,time,mag,merr,tbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nbins,bin,nmax,j
      parameter(nmax=650000)
      real time(npt),mag(npt),merr(npt),tbin,tmin,tmax,ltime,
     .     avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .     var(nmax),ep(nmax),s,p,avgs(nmax),sigcut


      sigcut=3.0

C     assume time is in days and tbin is in minutes

      do 5 i=1,nmax
         avgm(i)=0.0
         avge(i)=0.0
         avgt(i)=0.0
         abin(i)=0.0
 5    continue

      tmin= 99.9e30
      tmax=-99.9e30
      do 10 i=1,npt
         tmin=min(tmin,time(i))
         tmax=max(tmax,time(i))
 10   continue

      ltime=tmax-tmin
      nbins=int(ltime*24.0*60.0/tbin+0.5)+1

      do 20 i=1,npt
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
         if(bin.gt.nmax) write(6,*) "WARNING, nmax too small in bindt"
         avgm(bin)=avgm(bin)+mag(i)/merr(i)
         avgt(bin)=avgt(bin)+time(i)/merr(i)
         avge(bin)=avge(bin)+1.0
         abin(bin)=abin(bin)+1.0/merr(i)
 20   continue

      do 21 i=1,nbins
         avgs(i)=avgm(i)/abin(i)
         ep(i)=0.
         var(i)=0.
         abin(i)=0
 21   continue

      do 22 i=1,npt
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
         s=mag(i)-avgs(bin)
         ep(bin)=ep(bin)+s
         p=s*s
         var(bin)=var(bin)+p
         abin(bin)=abin(bin)+1.0
 22   continue

      do 23 i=1,nbins
         var(i)=(var(i)-ep(i)**2/abin(i))/(abin(i)-1)
         stdev(i)=sqrt(var(i))
         if(abin(i).eq.0) stdev(i)=0.
c         write(6,*) avgs(i),stdev(i)
 23   continue

      do 6 i=1,nmax
         avgm(i)=0.0
         avge(i)=0.0
         avgt(i)=0.0
         abin(i)=0.0
 6    continue

      do 24 i=1,npt
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
         if(abs(mag(i)-avgs(bin)).lt.sigcut*stdev(bin)) then
            avgm(bin)=avgm(bin)+mag(i)/merr(i)
            avgt(bin)=avgt(bin)+time(i)/merr(i)
            avge(bin)=avge(bin)+1.0
            abin(bin)=abin(bin)+1.0/merr(i)
         endif
 24   continue

      j=0
      do 30 i=1,nbins
         if(abin(i).gt.0.0) then
            j=j+1
            avgm(i)=avgm(i)/abin(i)
            avgt(i)=avgt(i)/abin(i)
            avge(i)=(avge(i)**0.5)/abin(i)
            time(j)=avgt(i)
            mag(j)=avgm(i)
            merr(j)=avge(i)
c            write(6,*) time(j),mag(j),merr(j)
         endif
 30   continue
      npt=j
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bind(npt,time,mag,merr,bins)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,nmax,i,j,nc
      parameter(nmax=650000)
      real mag(nmax),magn(nmax),time(nmax),merr(nmax),timen(nmax),
     .  merrn(nmax),avgm,avgt,avge,we

      j=0
      avgm=0.
      avgt=0.
      avge=0.
      nc=0
      we=0.
      do 10 i=1,npt
        j=j+1
        avgm=avgm+mag(i)/merr(i)
        avgt=avgt+time(i)/merr(i)
        avge=avge+merr(i)
        we=we+1.0/merr(i)   
        if(j.eq.bins) then
          avgm=avgm/we
          avgt=avgt/we
          avge=avge/real(bins)
          nc=nc+1
          magn(nc)=avgm
          timen(nc)=avgt
          merrn(nc)=avge       
          j=0
          avgm=0
          avgt=0
          avge=0
          we=0.
        endif
 10   continue
  
      do 20 i=1,nc
        time(i)=timen(i)
        mag(i)=magn(i)
        merr(i)=merrn(i)
c        write(6,*) time(i),mag(i),merr(i)
  20  continue
      npt=nc

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

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotdata(npt,time,mag,merr,x1,x2,panx,pany,MOSTtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will plot the data (ie. time vs mag)
C     integer npt: Total number of data points
C     real time(npt): The x-axis data points
C     real mag(npt): The y-axis data points
C     real merr(npt): Error associated with mag()
C     real x1,x2: range of data to plot
C     panx,pany: PGPLOT panel to plot data in

      implicit none
      integer npt,panx,pany,i,sym
      real time(npt),mag(npt),merr(npt),x1,x2,mmin,mmax,tmin,tmax
      double precision MOSTtime
      character*80 xlabel

      sym=17
      if(npt.gt.1000) sym=-1
     

      call pgpanl(panx,pany)

      mmin= 99.9e30
      mmax=-99.9e30
      tmin= 99.9e30
      tmax=-99.9e30

      do 10 i=1,npt
         if(mag(i).gt.mmax) mmax=mag(i)
         if(mag(i).lt.mmin) mmin=mag(i)
         if(time(i).gt.tmax) tmax=time(i)
         if(time(i).lt.tmin) tmin=time(i)
 10   continue

      if(x1.ne.-1.0) then
         tmin=x1
      endif
      if(x2.ne.-1.0) then
         tmax=x2
      endif

c      mmax=0.022
c      mmin=-0.008
c      mmax=22.0
c      mmin=-5.0
      write(6,*) "mmin, mmax",mmin,mmax
      call pgwindow(tmin,tmax,mmax,mmin)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      write(xlabel,500) "HJD-",MOSTtime
 500  format(A4,F16.8)
c      call pglabel(xlabel,"mmag","")
      call pgptxt((tmin+tmax)/2.0,mmax+0.27*(mmax-mmin),0.0,0.5,xlabel)
      call pgptxt(tmin-0.05*(tmax-tmin),(mmax+mmin)/2,90.0,0.5,"mmag")
      call pgsch(2.0)
      call pgpt(npt,time,mag,sym)
      call pgsch(2.9)
c      call pgerrb(6,npt,time,mag,merr,1.0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function nyquest(npt,time)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ndt,i
      parameter(nmax=650000)
      real time(npt),dt(nmax),mindt,maxdt,mode,modecal,mean,std,stdev,
     .     sigcut

      sigcut=1.0

      ndt=0
      mindt= 99.9e30
      maxdt=-99.9e30
      mean=0.
      do 10 i=2,npt
         ndt=ndt+1
         dt(ndt)=time(i)-time(i-1)
         mindt=min(dt(ndt),mindt)
         maxdt=max(dt(ndt),maxdt)
         mean=mean+dt(ndt)
 10   continue
      mean=mean/real(ndt)
      std=stdev(ndt,dt,mean)
      mindt=max(mindt,mean-sigcut*std)
      maxdt=min(maxdt,mean+sigcut*std)

      mode=modecal(ndt,dt,mindt,maxdt)
      
c      write(6,*) "mean,mode,mindt,maxdt:",mean,mode,mindt,maxdt

      nyquest=1.0/(2.0*mode)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function modecal(npt,pts,dmin,dmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbins,i,j
      parameter(nbins=200)
      integer bins(nbins),maxbin,maxi
      real pts(npt),dmin,dmax,bspace,mode,x(nbins),y(nbins),minx,maxx,
     .     miny,maxy

      do 10 i=1,nbins
         bins(i)=0
 10   continue
      maxbin=-99

      modecal=300.0

      do 20 j=1,npt
C     calculating the bspacing between each bin
         bspace=(dmax-dmin)/real(nbins)
C     calculating the bin number
         i=(pts(j)-dmin)/bspace
         if((i.gt.0).and.(i.le.nbins)) then
c            write(6,*) i
            bins(i)=bins(i)+1
            if(bins(i).gt.maxbin) then 
               maxbin=bins(i)
               maxi=i
               mode=real(i)*bspace+dmin
            endif
         endif
 20   continue

      modecal=mode

c      write(6,*) "mode,nbin",modecal,maxbin
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function integrate(npt,time,mag,period,k,sc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function evaluates the fourier co-efficents though integral
C     approximations.  These parameters are then used a first guess for 
C     least square minimization
C
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     period - fixed period 
C     k - order of fourier series
C     sc - sine (=1) or cosine (=2) term
      implicit none
      integer nmax
      parameter (nmax=650000)
      integer npt,iord(nmax),i,k,sc
      real time(npt),mag(npt),period,fint,phase(nmax),pi

      pi=3.141592654

      call phasept(npt,time,mag,phase,period)
      call rqsort(npt,phase,iord)
      
      do 5 i=1,npt
         phase(i)=phase(i)-0.5
 5    continue
      
      if(k.lt.1) then
         fint=(phase(iord(1))+0.5)*mag(iord(1))
         do 10 i=2,npt
            fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))
 10      continue
      else
         if(sc.eq.2) then
            fint=(phase(iord(1))+0.5)*mag(iord(1))*cos(2.0*real(k)*pi*
     .           phase(iord(1)))
            do 20 i=2,npt
               fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))*
     .              cos(2.0*real(k)*pi*phase(iord(i)))
 20         continue
            fint=2.0*fint
         elseif(sc.eq.1) then
            fint=(phase(iord(1))+0.5)*phase(iord(1))*sin(2.0*real(k)*pi*
     .           mag(iord(1)))
            do 30 i=2,npt
               fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))*
     .              sin(2.0*real(k)*pi*phase(iord(i)))
 30         continue
            fint=2.0*fint
         else
            fint=-1.0
         endif
      endif

      integrate=fint

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine shaperm(npt,nfit,time,mag,merr,res,pers,avg,a,panx,
     .     pany,plot,w,nstar,MOSTtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Given fourier coefficents, this shape is removed from the data.
      implicit none
      integer npt,nfit,na,nstarmax,panx,pany,plot
      parameter(nstarmax=300)
      real time(npt),mag(npt),res(npt),pers(nstarmax),avg,a(nfit),
     .     merr(npt)
      double precision MOSTtime
      integer i,j,k,nstar,i2,cc
      real pi,w(nstarmax),temp,arg

      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar

      j=(na-nstar)/nstar
      do 10 i2=1,npt
         temp=avg
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
               arg=real((i-cc+2)/2)*a(cc-1)*time(i2)+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         res(i2)=mag(i2)-temp
 10   continue

      if(plot.eq.1) call plotdata(npt,time,res,merr,-1.0,-1.0,panx,pany,
     . MOSTtime)
      
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowfn(npt,nfit,time,mag,merr,period,per1,per2,avg,a,
     .   nb,nc,steps,panx,pany,nstar,pers)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Attempts to model the window function for the data by using shape
C     parameters from COSINEFIT.  A noiseless light curve is sampled 
C     using the TIME(npt) data.  The model data is then run through 
C     PDM to produce a periodogram that approximates the window 
C     function.
C     
      implicit none
      integer nmax,npt,nfmax,nfit,i,j,nb,nc,steps,panx,pany,nstarmax,
     .   nstar,i2,k,cc,na
      parameter(nmax=650000,nstarmax=300)
      real time(npt),mag(npt),per1,per2,avg,a(nfit),temp,res(nmax),
     .     bper,period,w(nstarmax),pi,btheta,merr(npt),pers(nstarmax),
     .     arg,sn,snlimit

      snlimit=3.6 !S/N for jmfourw (useless for window function)
      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar
      
      j=(na-nstar)/nstar
      
      write(6,*) "na,j",na,j
      do 10 i2=1,npt
         temp=0.0
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
               arg=real((i-cc+2)/2)*a(cc-1)*time(i2)+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         res(i2)=temp
 10   continue


      call jmfourw(npt,time,res,merr,per1,per2,steps,bper,btheta,sn,
     .   panx,pany,snlimit)
c      call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,panx,pany)
C      call plotph(npt,time,res,period,1,4)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(npt,time,mag,period,panx,pany)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     period - fixed period of the data
C     panx,pany - pgplot panel to plot light curve
      implicit none
      integer nmax
      parameter(nmax=650000)

      integer npt,panx,pany,i
      real time(npt),mag(npt),period,phase(nmax),mmax,mmin,ax1,ax2
      real temp,stdev,std,sigcut,mean

      mean=0.
      do 4 i=1,npt
         mean=mean+mag(i)
 4    continue
      mean=mean/real(npt)

      sigcut=5.0
      mmin= 99.9e30
      mmax=-99.9e30

      mmax=-10.0e10
      mmin= 10.0e10

      std=stdev(npt,mag,mean)
      call pgpanl(panx,pany)
      call phasept(npt,time,mag,phase,period)
      do 10 i=1,npt
         if((mag(i).gt.mmax).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
         if((mag(i).lt.mmin).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
 10   continue

      ax2=mmin-0.10*(mmax-mmin)
      ax1=mmax+0.10*(mmax-mmin)
      call pgwindow(0.0,1.0,ax1,ax2)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Phase","magnitude"," ")
	  call pgptxt(0.5,ax1+0.27*(ax1-ax2),0.0,0.5,"Phase")
      call pgptxt(-0.05,(ax1+ax2)/2,90.0,0.5,"mmag")
      if (npt.lt.1000) then
         call pgpt(npt,phase,mag,17)
      else
         call pgpt(npt,phase,mag,-1)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grelor(npt,time,mag,merr,m,per1,per2,steps,panx,pany)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine is a work in progress...
C
      implicit none
      integer npt,m,steps,maxm,funcN,maxstep,panx,pany,nmax
      parameter(maxm=20,maxstep=100000,nmax=650000)
      real time(npt),mag(npt),merr(npt),per1,per2,w1,w2,pi,px(maxstep),
     .     py(maxstep),ymax,ymin,tempr,zpt,fnorm,flux(nmax),pxl(2),
     .     pyl(2)
      double precision Wmult,Wm,w,temp1,nchoosek,sj(maxm),S,Sbfac,v,
     .     phshft,func,temp,ans,Om1,Om1w(maxstep),dsteps
      character*79 tempc
      integer nbin(maxm),i
      
      call pgpanl(panx,pany)
      call pgwindow(per1,per2,-8.0,15.0)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Period","Om1"," ")

      zpt=30.0
      pi=3.141592654
      pxl(1)=0.0
      pyl(1)=0.0

      w1=2.0*pi/per2
      w2=2.0*pi/per1

      tempr=0.
      do 5 i=1,npt
         flux(i)=10.0**((mag(i)-zpt)/-2.5)
         tempr=flux(i)+tempr
 5    continue
      fnorm=tempr/real(npt)
      do 6 i=1,npt
         flux(i)=flux(i)/fnorm
 6    continue

      temp1=nchoosek(npt+m-1,npt)
      v=dble(m-1)
      Om1=0.0D0
      dsteps=dble(steps)

      do 10 i=1,steps
         w=w2-(w2-w1)/dsteps*dble(i-1)
         call qsimp(funcN,0.0D0,1.0D0,ans,npt,time,flux,m,nbin,sj,w)
         ans=ans/(2.0D0*dble(pi)*v*temp1)
         temp=(w*dble(log(w2/w1)))**-1.0D0
         ans=ans*temp
         px(i)=real(2.0*pi/real(w))
         py(i)=real(log10(ans))
         write(tempc,*) px(i),py(i),"           "
         call ovrwrt(tempc,2)
         pxl(2)=px(i)
         pyl(2)=py(i)
         call pgline(2,pxl,pyl)
         pxl(1)=px(i)
         pyl(1)=py(i)
         Om1w(i)=ans
         Om1=Om1+ans
 10   continue
      
      Om1=Om1*(w2-w1)
      write(6,*) Om1, "Prob: ",Om1/(1+Om1)

      ymin= 99.9e30
      ymax=-99.9e30
      do 20 i=1,steps
         if(py(i).lt.ymin)ymin=Om1w(i)
         if(py(i).gt.ymax)ymax=Om1w(i)
 20   continue

c      tempr=(ymax-ymin)*0.2
c      call pgpanl(panx,pany)
c      call pgwindow(per1,per2,ymin-tempr,ymax+tempr)
c      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Period","Om1"," ")
c      call pgpt(steps,px,py,17)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function func(npt,time,flux,m,nbin,sj,w,phshft)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,m,nbin(m)
      real time(npt),flux(npt),fnorm
      double precision sj(m),S,Wm,w,phshft,Sbfac,Wmult

      call bindata(npt,time,flux,w,m,nbin,sj,phshft)
      S =Sbfac(m,nbin,sj)
      Wm=Wmult(npt,m,nbin,1)

      func=S/Wm

      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Sbfac(m,nbin,sj)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer m,nbin(m),i
      double precision sj(m),ans

      ans=0
      do 10 i=1,m
         if(sj(i).gt.0.0d0) ans=ans+sj(i)**(-nbin(i))
 10   continue

      Sbfac=ans

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(npt,time,flux,w,m,nbin,sj,phshft)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,m,nbin(m),nmax,i,tmpi,maxm
      parameter(nmax=650000,maxm=20)
      real time(npt),flux(npt),phase(nmax),period,pi,rbin(maxm),zpt
      double precision w,sj(m),temp,phshft
      pi=3.141592654
      zpt=30.0

      period=2.0*pi/real(w)
      call phasept(npt,time,flux,phase,period)

      do 5 i=1,m
         rbin(i)=0.0
 5    continue

      temp=real(phshft)
      do 10 i=1,npt
         phase(i)=phase(i)+temp
         phase(i)=phase(i)-int(phase(i))
         tmpi=int(m*phase(i))+1
         rbin(tmpi)=rbin(tmpi)+flux(i)
 10   continue
      
      temp=dble(npt)/dble(m)

      do 20 i=1,m
         nbin(i)=int(rbin(i)+0.5)
         sj(i)=dble(nbin(i))/temp
 20   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function nchoosek(n,k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,k
      integer i,nbin(2)
      double precision Wmult

      nbin(1)=n-k
      nbin(2)=k
      
      nchoosek=Wmult(n,2,nbin,0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Wmult(npt,m,nbin,which)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,m,nbin(m),i,cont,maxm
      parameter(maxm=20)
      integer tbin(m),tpt,maxbin,tempmax,which
      double precision ans,lexpo,x,dbm

      do 10 i=1,m
         tbin(i)=nbin(i)
 10   continue
      tpt=npt
      
      cont=m+1
      if(which.eq.1) then
         dbm=dble(m)
         lexpo=dble(npt)
         x=1.0D2/log10(dbm)
         ans=1.0D0/(dbm**x)
         lexpo=lexpo-x
      else
         ans=1.0
      endif
      do while (cont.gt.0)
         tempmax=tbin(1)
         maxbin=1
         do 20 i=2,m
            if(tbin(i).gt.tempmax) then
               maxbin=i
               tempmax=tbin(i)
            endif
 20      continue
         if(tbin(maxbin).gt.1) then
            ans=ans*dble(tpt)/dble(tbin(maxbin))
            tpt=tpt-1
            tbin(maxbin)=tbin(maxbin)-1
            if((which.eq.1).and.(ans.gt.1.0d0).and.(lexpo.gt.x))then
               ans=ans/(dble(m)**x)
               lexpo=lexpo-x
            endif
         else
            cont=0
         endif
C         write(6,*) tbin,tpt,ans,lexpo
      enddo

      do 30 i=2,tpt
         ans=ans*dble(i)
         if((which.eq.1).and.(ans.gt.1.0d0).and.(lexpo.gt.x))then
            ans=ans/(dble(m)**x)
            lexpo=lexpo-x
         endif
 30   continue
      ans=ans/(dble(m)**lexpo)

      Wmult=ans
         
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pdm(npt,time,mag,nb,nc,per1,per2,steps,bper,btheta,
     .     panx,pany)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt -  number of data points
C     time - time of observations
C     mag  - observations
C     nb,nc - bin struction (classically=5,2)
C     per1,per2 - period range to search
C     steps - number of periods to search in range
C     panx,pany - which pgplot panel to plot the results in
      implicit none
      integer nmax,stepmax
      parameter(nmax=650000,stepmax=1000)

      integer npt,nb,nc,steps,panx,pany
      real time(npt),mag(npt),per1,per2,bper,period
      
      integer i,j
      real sumq,sumq2,sigma2,theta,phase(nmax),FTHETA,btheta,temp,
     .     nx(stepmax),ny1(stepmax),ny2(stepmax),px(3),py(3),ave,
     .     adev,sdev,skew,curt,perc,percold
      character*80 tmpc

      if(steps.lt.stepmax) then
         write(6,*) "Increasing number of steps to 2000"
         steps=stepmax
      endif

C     Sums for computing standard deviation
      sumq=0
      sumq2=0

      do 5 i=1,stepmax
         ny1(i)=-99.9e10
         ny2(i)= 99.9e10
 5    continue

C     calculate the variance
c      do 10 i=1,npt
c         sumq=sumq+mag(i)
c         sumq2=sumq2+mag(i)**2
c 10   continue
c      sigma2= (sumq2-sumq**2/real(npt))/real(npt-1)
      call moment(mag,npt,ave,adev,sdev,sigma2,skew,curt)
      
      percold=0.
      btheta=99.9e10
      temp=per2-per1
      if(temp.lt.0) temp=abs(temp)
      do 20 i=1,steps
         perc=100.0*real(i)/real(steps)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
            percold=perc
            call ovrwrt(tmpc,2)
         endif
 500     format(A13,F6.1)
         j=1+real(stepmax)/real(steps)*(i-1)
         period=per1+real(i)*(temp)/real(steps)
         period=1.0/period
         call phasept(npt,time,mag,phase,period)
         theta=FTHETA(npt,phase,mag,period,nb,nc,sigma2)
         if(theta.eq.-1) theta=1.0
         nx(j)=1.0/period
         if(theta.gt.ny1(j)) ny1(j)=theta
         if(theta.lt.ny2(j)) ny2(j)=theta
         if(theta.lt.btheta) then
            btheta=theta
            bper=1.0/period
         endif
 20   continue
c      write(6,*) bper,btheta
c      call pgpage()
      call pgpanl(panx,pany)
      call pgwindow(per1,per2,btheta,1.00)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("Period","Theta"," ")
      call pgbbuf()
      do 30 i=2,stepmax
         px(1)=nx(i-1)
         py(1)=ny2(i-1)
         px(2)=nx(i)
         px(3)=nx(i)
         py(2)=ny1(i)
         py(3)=ny2(i)
         call pgline(3,px,py)
 30   continue
      call pgebuf()

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
      real function FTHETA(npt,phase,mag,period,nb,nc,sigma2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     return theta statistic
C     npt - number of data points
C     phase - phased data points
C     mag - the observations
C     period - fixed period for data
C     nb,nc - bin structure
C     sigma2 - variance of data
      implicit none
      integer npt,nb,nc
      real phase(npt),mag(npt),period,sigma2,temp,fsvar

      temp=FSVAR(npt,phase,mag,period,nb,nc)

      if(temp.eq.-1.0) then
         FTHETA=-1.0
      else
         FTHETA=temp/sigma2
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function FSVAR(npt,phase,mag,period,nb,nc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nb,nc
      real phase(npt),mag(npt),period

      real temp1,temp2,temp3,sum(0:10,0:10),sum2(0:10,0:10),var
      integer nbr,ncr,shift,group,ptpbin(0:10,0:10),j

      temp1=0
      temp2=0
      nbr=nb
      ncr=nc

      do 10 shift=0,nc-1
         do 5 group=0,nb-1
            ptpbin(shift,group)=0
            sum(shift,group)=0
            sum2(shift,group)=0
 5       continue
 10   continue
      do 30 shift=0,nc-1
         do 20 j=1,npt
            temp3=phase(j)-shift/nbr/ncr+1
            group=int((temp3-int(temp3))*nbr)
            sum(shift,group)=sum(shift,group)+mag(j)
            sum2(shift,group)=sum2(shift,group)+mag(j)**2
            ptpbin(shift,group)=ptpbin(shift,group)+1
 20      continue
 30   continue
      do 50 shift=0,nc-1
         do 40 group=0,nb-1
            if (ptpbin(shift,group).le.1) goto 60
            var=(sum2(shift,group)-sum(shift,group)**2/real(
     .           ptpbin(shift,group)))/real(ptpbin(shift,group)-1)
            temp1=temp1+(ptpbin(shift,group)-1)*var
            temp2=temp2+ptpbin(shift,group)
 40      continue
 50   continue
      goto 70
 60   fsvar=-1
      goto 80
 70   fsvar=temp1/(temp2-nbr*ncr)
 80   return
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


ccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine updatecosinefit(a,rmavg,ma,nstar,aerr)
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
      parameter(nmax=650000,nfmax=2000,nfitm=2000,nstarmax=100)
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

ccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cosinefit(pers,res,id,xcoo,ycoo,avg,a,btheta,rmavg,
     .     plot,ma,aerr,panx,pany,w,nstar,fixfreq)
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
      integer nmax,nfmax,nfit,nfitm,b,nstar,nstarmax,plot,endfit,panx,
     .     pany,fixfreq
      parameter(nmax=650000,nfmax=2000,nfitm=2000,nstarmax=300)
      real w(nstarmax), x(nmax),y(nmax),sig(nmax),a(nfmax),chisq,tx,ty,
     .     ph(nmax),alpha(nfmax,nfmax),covar(nfmax,nfmax),alamda,avg,
     .     avgc,per,ty2(nmax),tx2(nmax),avgbin(10),max,min,amp,bad,
     .     xcoo,ycoo,xi(nfitm*2),ctrl(12),f,pi,maxamp,
     .     fit1a,temp(2),integrate,btheta,rmavg,pers(nstarmax),arg,
     .     twopi,res(nmax),ochi,alamhi,chitol,aerr(nfmax)
      integer(kind=1) now(3)
      integer i,n,ia(nfmax),ma,nca, j, k,avgn(10),id,seed,status,cc,j1,
     .     i2,niter
      character*50 filename,starname
      logical loop
      common /data/ x,y,sig,n,nfit
c      external fit1a

      alamhi=100000.
      chitol=10e-6

      maxamp=2.0
      pi=3.141592654
      twopi=2*pi
      do 4 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
c         write(6,*) "frequency",w(i),pers(i)
 4    continue

      loop=.true.
      i=0
      avg=0.
      max=-1.0e25
      min=1.0e25
      do 34 i=1,n
         if(y(i).gt.max)max=y(i)
         if(y(i).lt.min)min=y(i)
c         ph(i)=x(i)/pers(1) - int(x(i)/pers(1))
 34   continue
      
      avg=integrate(n,x,y,pers(1),0,1)      
      
      do 5 i=1,n
         y(i)=y(i)-avg
 5    continue
      
      ma=nstar*nfit*2+nstar
      do 10 i=1,ma
         ia(i)=1.0
 10   continue
      nca=nfmax

      alamda=-1.0

      j=(ma-nstar)/nstar
      do 30 k=nstar,nstar
         cc=(j+1)*(k-1)+2
         a(cc-1)=w(k)
         per=pers(k)
         do 31 i=cc,cc+j-2,2
            if(nstar.gt.1) then
               temp(1)=integrate(n,x,res,per,(i-cc+2)/2,1)
               temp(2)=integrate(n,x,res,per,(i-cc+2)/2,2)
            else
               temp(1)=integrate(n,x,y,per,(i-cc+2)/2,1)
               temp(2)=integrate(n,x,y,per,(i-cc+2)/2,2)
            endif
            a(i+1)=atan(-temp(1)/temp(2))
            a(i)=temp(2)/cos(a(i+1))
            a(i+1)=a(i+1)+pi
 31      continue
 30   continue

      do 32 k=1,nstar
         cc=(j+1)*(k-1)+2
C     THIS LINE HANDLES FIXED PERIODS
         if(fixfreq.gt.0) ia(cc-1)=0
         do 33 i=cc,cc+j-2,2
            if(a(i).lt.0.0) then
               a(i)=abs(a(i))
               a(i+1)=a(i+1)+pi
            endif
            if(a(i+1).lt.0) a(i+1)=2*pi+a(i+1)
            temp(1)=a(i+1)/(2.0*pi)
            b=int(temp(1))
            a(i+1)=a(i+1)-real(b)*2.0*pi
 33      continue
 32   continue

c      goto 50

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
C         write(6,*) avgc

C     look for negitive amplitudes!
CC         j=(ma-nstar)/nstar
CC         do 38 k=1,nstar
CC            cc=(j+1)*(k-1)+2
CC            do 39 i=cc,cc+j-2,2
CC               if(a(i).lt.0.0) then
CC                  a(i)=abs(a(i))
CC                  a(i+1)=a(i+1)+pi
CC               endif
CC               if (a(i+1).lt.0) a(i+1)=twopi+a(j+1)
CC               b = int(a(i+1)/twopi)
CC               a(i+1)=a(i+1)-b*twopi
CC 39         continue
CC 38      continue
         if(chisq.gt.ochi) then
            endfit=0
         elseif(abs(chisq-ochi).lt.chitol) then
            endfit=endfit+1
         endif
      enddo
c 20   continue
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

      if(plot.eq.1) then
         if(nstar.eq.1) then 
            call plotph(n,x,y,per,panx,pany)
         else
            call plotph(n,x,res,per,panx,pany)
         endif
         call pgsci(2)
         call pgpt(1000,tx2,ty2,-1)
         call pgsci(1)
      endif
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
         write(6,501) "Fre", a(cc-1)/twopi,aerr(cc-1)
         do 43 i=cc,cc+j-2,2
            aerr(i)=sqrt(covar(i,i))
            aerr(i+1)=sqrt(covar(i+1,i+1))
 43      continue
         write(6,500) "Amp", (a(i)  ,i=cc,cc+j-2,2)
         write(6,500) "err", (sqrt(abs(covar(i,i))),    i=cc,cc+j-2,2)
         write(6,500) "Psi", (a(i+1),i=cc,cc+j-2,2)
         write(6,500) "err", (sqrt(abs(covar(i+1,i+1))),i=cc,cc+j-2,2)
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
      parameter(nstarmax=300)
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

C      include "pikaia.f"

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE qsimp(funcN,a,b,s,npt,time,flux,m,nbin,sj,w)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Adapted from Numerical Recipes in Fortran 77 Second Edition
C
C     This routine uses the extended trapezoidal rule for numerically
C     evalutating a definite integral and defined by the function
C     FUNC(X,N)
C
      INTEGER JMAX,npt,m,nbin(m)
      INTEGER FUNCN
      DOUBLE PRECISION a,b,func,s,EPS,sj(m),w
      real time(npt),flux(npt)
      PARAMETER (EPS=1.e-2, JMAX=20)
C     USES trapzd
      INTEGER j
      DOUBLE PRECISION os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX
         call trapzd(funcN,a,b,st,j,npt,time,flux,m,nbin,sj,w)
         s=(4.*st-ost)/3.
         if (j.gt.5) then
            if (abs(s-os).lt.EPS*abs(os).or.
     *           (s.eq.0..and.os.eq.0.)) return
         endif
         os=s
         ost=st
c         write(6,*) s
 11   enddo
      pause 'too many steps in qsimp'
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE trapzd(funcN,a,b,s,n,npt,time,flux,m,nbin,sj,w)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Adapted from Numerical Recipes in Fortran 77 Second Edition
C
C     Implements the extended trapezoidal rule for evaluating a definite
C     Integral.
C
      INTEGER n,FUNCN,npt,m,nbin(m)
      DOUBLE PRECISION a,b,s,func,sj(m),w
      real time(npt),flux(npt)
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
         s=0.5*(b-a)*(func(npt,time,flux,m,nbin,sj,w,a)+
     .        func(npt,time,flux,m,nbin,sj,w,b))
      else
         it=2**(n-2)
         tnm=it
         del=(b-a)/tnm
         x=a+0.5*del
         sum=0.
         do 11 j=1,it
            sum=sum+func(npt,time,flux,m,nbin,sj,w,x)
            x=x+del
 11      enddo
         s=0.5*(s+(b-a)*sum/tnm)
      endif
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
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (6,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (6,2) output, char(13), char(13)
            write (6,2) output, char(13), char(13)
            write (6,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (6,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (6,3) line
    3    format (a)
      else
         write (6,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (6,2) line, char(13), char(13)
      end if
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,ia(ma),npc,ndat,MMAX 
      REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat) 
      EXTERNAL funcs2 
      PARAMETER (MMAX=50)

      INTEGER i,j,k,l,m,mfit 
      REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX) 
      mfit=0 
      do 11 j=1,ma 
         if(ia(j).ne.0) mfit=mfit+1 
 11   continue
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted'
      do 13 j=1,mfit 
         do 12 k=1,mfit 
            covar(j,k)=0. 
 12      continue
         beta(j)=0. 
 13   continue
      do 17 i=1,ndat
         call funcs2(x(i),afunc,ma) 
         ym=y(i) 
         if(mfit.lt.ma) then 
            do 14 j=1,ma 
               if(ia(j).eq.0) ym=ym-a(j)*afunc(j) 
 14         continue
         endif 
         sig2i=1./sig(i)**2 
         j=0 
         do 16 l=1,ma 
            if (ia(l).ne.0) then 
               j=j+1 
               wt=afunc(l)*sig2i 
               k=0 
               do 15 m=1,l 
                  if (ia(m).ne.0) then 
                     k=k+1 
                     covar(j,k)=covar(j,k)+wt*afunc(m) 
                  endif 
 15            continue
               beta(j)=beta(j)+ym*wt 
            endif 
 16      continue 
 17   continue 
      do 19 j=2,mfit 
         do 18 k=1,j-1 
            covar(k,j)=covar(j,k) 
 18      continue
 19   continue 
      call gaussj(covar,mfit,npc,beta,1,1) 
      j=0 
      do 21 l=1,ma 
         if(ia(l).ne.0) then 
            j=j+1 
            a(l)=beta(j) 
         endif 
 21   continue
      chisq=0.
      do 23 i=1,ndat 
         call funcs2(x(i),afunc,ma) 
         sum=0. 
         do 22 j=1,ma 
            sum=sum+a(j)*afunc(j) 
 22      continue 
         chisq=chisq+((y(i)-sum)/sig(i))**2 
 23   continue 
      call covsrt(covar,npc,ma,ia,mfit)
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE spline(x,y,n,yp1,ypn,y2) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER n,NMAX 
      REAL yp1,ypn,x(n),y(n),y2(n) 
      PARAMETER (nmax=650000) 
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

c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end
