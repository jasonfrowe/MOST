      program crosstalkfixer
      implicit none
      integer nmax,nunit,npt
      parameter(nmax=400000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),stdsky(nmax),
     .   maglow,maghi,gain,rdnoise,skydif(nmax),rmavg
      double precision dtime(nmax)
      character*80 filename
      
      gain=6.1 !e-/ADU
      rdnoise=1.1 !ADU
c      filename="hd209458_2004_01fs.dat"
c      filename="junk02.dat"
C      filename="test_01fs.dat"
      filename="hd209458_ap_01fcsacutd.dat"
      write(6,*) "filename:",filename
      
      nunit=10
      open(unit=nunit,file=filename,status='old',err=901)
      
C     hard cutoffs when reading in data
      maglow=-10.0
      maghi =13.0

C     call subroutine to read in data. 
      call readdataskyfit(nunit,npt,time,dtime,mag,merr,
     .     sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,
     .     stdsky,maglow,maghi)
      close(nunit)
      
      call pgopen('?') !query for plotting device
      call pgpage() !prepare for multiple pages
c      call windowsetup(xb1,xb2,yb1,yb2) !setup a square plot
c      call pgvport(xb1,xb2,yb1,yb2) !make a square plotting surface 
      call plotskys(npt,sky,stdsky,skydif,rdnoise,gain)
      call pgpage()
      call plotskydif(npt,time,mag,merr,skydif,stdsky,sky)
     
      rmavg=0.
      filename="test.dat"
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,stdsky,rmavg)
     
      call pgclos()
     
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotskydif(npt,time,mag,merr,skydif,stdsky,sky)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i,ii,bnpt,nlo,nhi,nphase,nfold,nptb
      parameter(nmax=400000)
      real time(npt),mag(npt),skydif(npt),ploty(nmax),stdsky(npt),ntest,
     .   tbin,bx(nmax),by(nmax),be(nmax),merr(npt),period,mindate,
     .   phase(nmax),bphase(nmax),bskydif(nmax),bmerr(nmax),bsky(nmax),
     .   sky(npt),bs(nmax),skydfind,magcor

      write(6,*) "calling rphase"
      period=0.0704256
      call rphase(npt,time,phase,period,nphase,mindate)
      nfold=50
C     if you want to do any binning, do it here.
      nptb=0
      do 4 i=1,npt
         ii=ntest(bskydif(i))
         if((ii.eq.0))then!.and.(sky(i).lt.4000.0))then
            nptb=nptb+1
            bphase(nptb)=phase(i)
            bskydif(nptb)=skydif(i)
            bmerr(nptb)=merr(i)
            bsky(nptb)=sky(i)
         endif
 4    continue
      call foldnremove(npt,time,skydif,nptb,bphase,bskydif,bmerr,
     .   nfold,mindate,period)

      bnpt=0
      do 5 i=1,npt
         ii=ntest(skydif(i))
         if(ii.eq.0)then
            bnpt=bnpt+1
            bx(bnpt)=time(i)
            by(bnpt)=skydif(i)
            be(bnpt)=merr(i)
            bs(bnpt)=sky(i)
         endif
 5    continue
      write(6,*) "testing",bnpt,npt
 
      tbin=101.!5.0
      nlo=3
      nhi=3
      write(6,*) "binning data"
      call newbindt(bnpt,bx,by,be,bs,tbin,nlo,nhi)
      
      open(unit=12,file="skydif.dat")
      do 6 i=1,bnpt
         write(12,*) bx(i),by(i)
 6    continue
      close(12)
      
      write(6,*) "getting plot data"
      do 10 i=1,npt
         ploty(i)=2500.0*mag(i)
c         ii=0
c         ii=ntest(skydif(i))
c         if((stdsky(i).lt.15).and.(ii.eq.0))write(6,*) mag(i),skydif(i)
         magcor=skydfind(time(i),mag(i),bnpt,bx,by)
         mag(i)=mag(i)-skydif(i)/2500.0
 10   continue

      write(6,*) "plotting"
      call pgwindow(2060.,2065.,100.,-100.)!set window scale
C      call pgwindow(2.,10.,100.,-100.)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)!boxes and ticks
      call pglabel("JD","\gD\gs","") !label axes
      call pgsci(3) 
      call pgpt(npt,time,ploty,17)
      call pgsci(2)
c      call pgpt(bnpt,bx,by,17)
      call pgline(bnpt,bx,by)
      call pgsci(1)

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function skydfind(time,mag,bnpt,bx,by)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer bnpt,i,j,nm,np
      real time,mag,bx(bnpt),by(bnpt),mdif,pdif,dif,m,b
      
      nm=1
      np=1
      mdif=bx(1)-time
      pdif=bx(1)-time
      do 11 j=2,bnpt
         dif=bx(j)-time
c         write(6,*) dif,mdif,pdif,time
         if(dif.le.0)then
            if(abs(dif).lt.abs(mdif))then
               mdif=dif
               nm=j
            endif
         else
            if(abs(dif).lt.abs(pdif))then
               pdif=dif
               np=j
            endif
         endif
 11   continue
c      write(6,*) dif,mdif,pdif,time,bx(1),nm,np
c      read(5,*)
      m=(by(np)-by(nm))/(bx(np)-bx(nm))
      b=by(np)-m*bx(np)

      skydfind=m*time+b
      if(bx(np)-bx(nm).eq.0.) skydfind=by(np)
c      write(6,*) time,mag,skydfind,bx(np),by(np)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine newbindt(npt,time,mag,merr,sky,tbin,nlow,nhigh)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
c     bins data into equal time bins
      implicit none
      integer npt,i,nbins,bin,nmax,j,ntemp,ntemp2,nlow,nhigh
      parameter(nmax=400000)
      real time(npt),mag(npt),merr(npt),tbin,tmin,tmax,ltime,
     .     avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .     var(nmax),ep(nmax),s,p,avgs(nmax),sigcut,xt(nmax),yt(nmax),
     .     zt(nmax),xt2(nmax),yt2(nmax),zt2(nmax),mode,modecal,sky(npt),
     .     st(nmax)


      sigcut=1.0

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

C     reject high and low values
      ntemp2=0
      if((nlow.gt.0).or.(nhigh.gt.0)) then
C     this will probably be really slow.. ugh.
         do 30 i=1,nbins
            ntemp=0
            do 31 j=1,npt
C     calculate bin number here.
               bin=int(real(nbins)*(time(j)-tmin)/ltime)+1
               if(bin.eq.i) then
                  ntemp=ntemp+1
                  xt(ntemp)=time(j)
                  yt(ntemp)=mag(j)
                  zt(ntemp)=merr(j)
                  st(ntemp)=sky(j)
               endif
 31         continue
c            if(ntemp.gt.0) then
c               mode=modecal(ntemp,st,0.0,800.)
c               write(6,*) "mode: ",ntemp,xt(1),mode
c            endif
C            read(5,*)
            call rejhilow(ntemp,xt,yt,zt,nlow,nhigh)
            do 32 j=1,ntemp
               ntemp2=ntemp2+1
               xt2(ntemp2)=xt(j)
               yt2(ntemp2)=yt(j)
               zt2(ntemp2)=zt(j)
 32         continue
 30      continue
         do 33 i=1,ntemp2
            time(i)=xt2(i)
            mag(i)=yt2(i)
            merr(i)=zt2(i)
 33      continue
         npt=ntemp2
      endif

      do 20 i=1,npt
         bin=int(real(nbins)*(time(i)-tmin)/ltime)+1
         if(bin.le.0) write(6,*) "wooaah... bin is stupid",bin
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

      do 6 i=1,npt
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
      do 40 i=1,nbins
         if(abin(i).gt.0.0) then
            j=j+1
            avgm(i)=avgm(i)/abin(i)
            avgt(i)=avgt(i)/abin(i)
            avge(i)=(avge(i)**0.5)/abin(i)
            time(j)=avgt(i)
            mag(j)=avgm(i)
            merr(j)=avge(i)
         endif
 40   continue
      npt=j
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine foldnremove(npt,time,mag,nptb,bphase,bmag,
     .   bmerr,nfold,mindate,period)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none 
      integer npt,nfold,i,j,nptb,nmax,nb,n
      parameter(nmax=400000)
      real time(npt),mag(npt),bphase(nptb),bmag(nptb),mindate,yp1,yp2,
     .   y2(nmax),bx(nmax),by(nmax),phase,pshift,y,period,perc,percold,
     .   ymin,ymax,avgx(nmax),avgy(nmax),avgz(nmax),bz(nmax),sig,
     .   bmerr(nmax),pbin
      character*80 tmpc,ans
      
      pbin=0.005
c 5    write(6,*) "Enter phase size for binning [0.005]"
c      read(5,501) ans
c 501  format(A80)
c      if(ans.eq." ") then
c         pbin=0.005
c      else
c         read(ans,*,err=5) pbin
c      endif
c      if((pbin.le.0.0).or.(pbin.ge.1.0))then
c         write(6,*) "That answer makes no sense"
c         goto 5
c      endif
             
      
      percold=0.
      
c      call pgopen('?')
      
      do 10 i=1,npt
c         call pgeras()
      
         perc=100.0*real(i)/real(npt)
         if(perc-percold.gt.0.1) then
            write(tmpc,500) "Percent done ", perc
 500        format(A13,F6.1)
            percold=perc
            call ovrwrt(tmpc,2)
         endif
c      write(6,*) "point number:",i,time(i),mag(i)
C     pivot point is currect dphase
C     we want all data +/- nfold/2
         nb=0
         phase=(time(i)-mindate)/period
         do 11 j=1,nptb
            if((bphase(j).gt.phase-real(nfold)/2).and.
     .         (bphase(j).lt.phase+real(nfold)/2)) then
               nb=nb+1
               bx(nb)=bphase(j)-int(bphase(j))
               by(nb)=bmag(j)
               bz(nb)=bmerr(j) 
            endif
 11      continue
C        now shift data so important phase is in the middle
         phase=phase-int(phase)
         pshift=phase-0.5
c         write(6,*) "pshift:",pshift
         do 12 j=1,nb
            bx(j)=bx(j)-pshift
            if(bx(j).lt.0.0) bx(j)=bx(j)+1.0
            if(bx(j).gt.1.0) bx(j)=bx(j)-1.0
 12      continue
C        sort data by phase
c         write(6,*) "Sorting.."
c         call sort2(nb,bx,by)
c         do 13 j=1,nb
c            write(6,*) j,bx(j),by(j)
c 13      continue
C        set up spline 
c         yp1=1.0e30  
c         yp2=1.0e30
c         call plotpoints(nb,bx,by,17,1,0.0,1.0,0.0,0.0)
         
cc         write(6,*) "Setting up spline"
c         call spline(bx,by,nb,yp1,yp2,y2)
cC        get interpolated point
cc         write(6,*) "Splining.."
c         call splint(bx,by,y2,nb,0.5,y)
         y=0
         n=0
         do 14 j=1,nb
            if(abs(0.5-bx(j)).lt.pbin) then
               y=y+by(j)
               n=n+1
               avgx(n)=bx(j)
               avgy(n)=by(j)
               avgz(n)=bz(j)
            endif
 14      continue
         
         sig=3.0
         call sigclip(n,avgx,avgy,avgz,sig)
         
         y=0
         do 15 j=1,n
            y=y+avgy(j)
 15      continue
         y=y/real(n)

c         call pgsci(3)
c         call pgpt1(0.5,y,10)
c         call pgsci(1)
         mag(i)=mag(i)-y
c         write(6,*) time(i),mag(i),phase,y
 10   continue
 
c      call pgclos()

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rphase(npt,time,phase,period,nphase,mindate)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nphase
      real period,time(npt),phase(npt),mindate
      
C     find the first date in the data set
      mindate=99.9d30
      do 5 i=1,npt
         mindate=min(mindate,time(i))
 5    continue
      
C     find phase of date (not integer chopping)

      nphase=-9999
      do 10 i=1,npt
         phase(i)=(time(i)-mindate)/period
c         write(6,*) time(i),phase(i)     
C        nphase is the larger integer phase found
         nphase=max(nphase,int(phase(i)))
 10   continue
 
      return
      end      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotskys(npt,sky,stdsky,skydif,rdnoise,gain)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nmax
      parameter(nmax=400000)
      real sky(npt),stdsky(npt),dxmin,dxmax,dymin,dymax,plotx(nmax),
     .   ploty(nmax),rdnoise,gain,lplot(nmax),x,y,skydif(npt)
 
      do 9 i=1,npt
         plotx(i)=rdnoise*gain+Sqrt(sky(i)*gain)
         x=plotx(i)
C     solution found using Mathematica.
c         y=-7.327812155995344 + 1.3665695990093576*x + 
c     .      0.00243344878334303*x**2
         y=11.615977709984893 + 0.37497344890404505*x + 
     .      0.0025687428240144337*x**2 + 8.3974781407879e-6*x**3
         ploty(i)=stdsky(i)*gain
         lplot(i)=y
         skydif(i)=ploty(i)-y
 9    continue
 
      dxmin=plotx(1)
      dxmax=plotx(1)
      dymin=ploty(1)
      dymax=ploty(1)
      do 10 i=2,npt
         dxmin=min(plotx(i),dxmin)
         dxmax=max(plotx(i),dxmax)
         dymin=min(ploty(i),dymin)
         dymax=max(ploty(i),dymax)
 10   continue
     
      call pgwindow(dxmin,dxmax,dymin,dymax)!set window scale
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)!boxes and ticks
      call pglabel("Expected \gs","Measured \gs","") !label axes 
      call pgpt(npt,plotx,ploty,-1)
      call pgsci(2)
      call pgpt(npt,plotx,lplot,1)
      call pgsci(1)
      
      return
      end
 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdataskyfit(nunit,npt,time,dtime,mag,merr,
     .     sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,
     .     stdsky,maglow,maghi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     id - id of star (returned)
C     xcoo,ycoo - co-ordinates of star (returned)
C     time,mag,merr - the data and associated error.
      implicit none
      integer nmax,nunit
      parameter(nmax=400000)
      integer npt,i,ntest,rn
      real time(nmax),mag(nmax),merr(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),maglow,maghi,stdsky(nmax)
      double precision dtime(nmax)

      i=1

   9  continue
  10  read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i),sky(i),
     .     xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),ftotal(i),
     .     aflux(i),nap(i),itime(i),stdsky(i)
        time(i)=real(dtime(i))
c         if(sky(i).gt.6000.0) goto 10
        merr(i)=abs(merr(i))
c  10  read(10,505,err=9,end=20) time(i),amag(i),mag(i),merr(i),sky(i),
c     .     xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i)
        if(merr(i).le.0.0) goto 10
c        if(fy(i).lt.0.5) goto 10
c        if(fx(i).lt.0.5) goto 10
        if((mag(i).lt.maglow).or.(mag(i).gt.maghi)) goto 10
C     Check for NaNs
        rn=ntest(mag(i))
        if(rn.eq.1) goto 10
        rn=ntest(sky(i))
        if(rn.eq.1) goto 10

        i=i+1
      goto 10
 20   continue

C     old data format.
c 505  format(F10.5,1X,3(F7.4,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
c     .     1X,F7.3,1X,F9.2)
C             13.8/11.6
 510  format(F13.8,1X,2(F9.6,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F9.2),1X,F8.3,1X,F6.2,1X,F8.4)

      npt=i-1
      write(6,*) "Points read: ",npt

      return
      end