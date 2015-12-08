C This code is property of Jason Rowe (2005)
C
C Contents
C
C avgrm
C sigclip
C binp
C binsig
C getspline
C modecal
C bind
C bindt
C rejhilow
C stdev
C sort
C rqsort




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine avgrm(npt,mag,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Removes average value from data points
C     This can help reduce numerical overflow for large data sets 
      implicit none
      integer npt
      real mag(npt),rmavg,avg

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
      subroutine sigclip(npt,time,mag,merr,sig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=600000,nitmax=50)
      real time(npt),mag(npt),merr(npt),sig,std,stdev,tmp1(nmax),
     .     tmp2(nmax),tmp3(nmax),mean,dum,omean

C     watch out for infinite loops
      niter=0

C     find mean of data set
      mean=0.
      do 4 i=1,npt
         mean=mean+mag(i)
 4    continue
      mean=mean/real(npt)
C     find standard dev. of data set.
      std=stdev(npt,mag,mean)

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
         else
            nclip=nclip+1
         endif
 10   continue
C     if nothing is clipped, no point in continuing
      if(nclip.eq.0) goto 15

C     save copy of old mean
      omean=mean

C     find mean of new data set.
      mean=0.
      do 5 i=1,ntmp
         mean=mean+tmp2(i)
 5    continue
      mean=mean/real(ntmp)
C     if mean does not change, we are done.
      if(mean.eq.omean) goto 15
C     find st. dev. of new data set
      std=stdev(ntmp,tmp2,mean)
C     now restart clipping...
      niter=niter+1
C     if we loop too much.. get out.
      if(niter.gt.nitmax) goto 15
      goto 6

 15   do 20 i=1,ntmp
         time(i)=tmp1(i)
         mag(i)=tmp2(i)
         merr(i)=tmp3(i)
 20   continue
      npt=ntmp

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine binp(npt,phase,mag,merr,bins)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,N,i,j
      parameter(N=600000)
      real phase(npt),tt(N),mag(npt),tn(N),merr(n),te(n)

      do 5 i=1,bins
         tn(i)=0.0
         tt(i)=0.0
         te(i)=0.0
 5    continue

      

      do 10 i=1,npt
         j=int(real(bins)*phase(i))+1
         if(j.le.bins)then
            tn(j)=tn(j)+1.0/merr(i)
            tt(j)=tt(j)+mag(i)/merr(i)
            te(j)=te(j)+1.0
         endif
 10   continue

      do 20 i=1,bins
         phase(i)=real(i)/real(bins)
         mag(i)=tt(i)/tn(i)
         merr(i)=(te(i)**0.5)/tn(i)
c         write(6,*) phase(i),mag(i),merr(i)
 20   continue
      npt=bins

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine binsig(npt,phase,mag,merr,sig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,N,i,nb,nmax,j
      parameter(N=1500,nmax=600000)
      integer iord(nmax),ntemp,ntemp2
      real phase(npt),mag(npt),merr(npt),sig,tt(n),tm(n),te(n),tn(n),
     .     xt(nmax),yt(nmax),zt(nmax),xt2(nmax),yt2(nmax),zt2(nmax)

      do 10 i=1,N
         tt(i)=0.
         tm(i)=0.
         te(i)=0.
         tn(i)=0
 10   continue

      call rqsort(npt,phase,iord)

      nb=1
      ntemp=0
      ntemp2=0
      do 40 i=1,npt
         tn(nb)=tn(nb)+1.0/merr(iord(i))
         te(nb)=te(nb)+1.0
         ntemp=ntemp+1
         xt(ntemp)=phase(iord(i))
         yt(ntemp)=mag(iord(i))
         zt(ntemp)=merr(iord(i))
          if(((te(nb)**0.5)/tn(nb).lt.sig).or.(i.eq.npt)) then
            call rejhilow(ntemp,xt,yt,zt,2,2)
            tn(nb)=0
            tt(nb)=0
            tm(nb)=0
            te(nb)=0
            do 50 j=1,ntemp
               tn(nb)=tn(nb)+1.0/zt(j)
               tt(nb)=tt(nb)+xt(j)/zt(j)
               tm(nb)=tm(nb)+yt(j)/zt(j)
               te(nb)=te(nb)+1.0
 50         continue
            nb=nb+1
            ntemp=0
         endif
 40   continue

      write(6,*) "nb:",nb

      do 30 i=1,nb
         phase(i)=tt(i)/tn(i)
         mag(i)=tm(i)/tn(i)
         merr(i)=(te(i)**0.5)/tn(i)
 30   continue
      npt=nb
      phase(nb-1)=(tt(nb)+tt(nb-1))/(tn(nb)+tn(nb-1))
      mag(nb-1)=(tm(nb)+tm(nb-1))/(tn(nb)+tn(nb-1))
      merr(nb-1)=((te(nb)+te(nb-1))**0.5)/(tn(nb)+tn(nb-1))
      npt=nb-1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getspline(npt,x,y,y2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
      integer npt
      real x(npt),y(npt),yp1,yp2,y2(npt)

      yp1=1.0e30
      yp2=1.0e30
      call spline(x,y,npt,yp1,yp2,y2)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function modecal(npt,pts,dmin,dmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
      implicit none
      integer npt,nbins,i,j
      parameter(nbins=50)
      integer bins(nbins),maxbin
      real pts(npt),dmin,dmax,bspace,mode

      do 10 i=1,nbins
         bins(i)=0
 10   continue
      maxbin=-99

      modecal=300.0

C     calculating the bspacing between each bin
      bspace=(dmax-dmin)/real(nbins)
      do 20 j=1,npt
C     calculating the bin number
         i=(pts(j)-dmin)/bspace
         if((i.gt.0).and.(i.le.nbins)) then
c            write(6,*) i
            bins(i)=bins(i)+1
            if(bins(i).gt.maxbin) then 
               maxbin=bins(i)
               mode=real(i)*bspace+dmin
            endif
         endif
 20   continue

      modecal=mode

c      write(6,*) "mode,nbin",modecal,maxbin
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bind(npt,x,y,yerr,bins,bsize,bx,by,berr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
C     bins data into bins of equal size
      implicit none
      integer npt,bins,bsize,nmax,i,j,ncut,k,nbins,bmax,niter,nitermax
      parameter(nmax=600000,bmax=3000)
      integer iord(nmax)
      real x(npt),y(npt),yerr(nmax),bx(bmax),by(bmax),mean,std,stdev,
     .     sigcut,omean,ty(nmax),meanx,stint,berr(bmax)

C     stop when mean converges
      stint=0.01
      
C     watch out for infinite loops in the sigma clipping
      nitermax=50

C     first we sort the data
      call rqsort(npt,x,iord)

      sigcut=2.0

      niter=0
      nbins=0
      do 10 i=1,npt,bsize
         nbins=nbins+1
         if(nbins.gt.bmax)write(6,*)"WARNING, increase BMAX"
         mean=0.
         k=0
         do 20 j=i,bsize+i-1
c            write(6,*) j,i,bsize,x(iord(j)),iord(j)
c            read(5,*)
            if(j.lt.npt) then
               k=k+1
               mean=mean+y(iord(j))
               ty(k)=y(iord(j))
            endif
 20      continue
         mean=mean/real(k)
 21      std=stdev(k,ty,mean)
c         write(6,*) "bind:",nbins,mean,std,k
         omean=mean
         k=0
         mean=0.
         meanx=0.
         ncut=0
         do 30 j=i,bsize+i-1
            if(j.lt.npt) then
               if(abs(y(iord(j))-omean).lt.sigcut*std) then
                  k=k+1
                  ty(k)=y(iord(j))
                  mean=mean+y(iord(j))
                  meanx=meanx+x(iord(j))
               else
                  ncut=ncut+1
               endif
            endif
 30      continue
         mean=mean/real(k)
         meanx=meanx/real(k)
         if((abs(mean-omean)/mean.lt.stint).or.(niter.gt.nitermax)) then
            bx(nbins)=meanx
            by(nbins)=mean
            berr(nbins)=std
c            write(6,*) nbins,bx(nbins),by(nbins),berr(nbins)
         else
c            write(6,*) "ni:",niter,mean,omean
            niter=niter+1
            goto 21
         endif
c         write(6,*) nbins,bx(nbins),by(nbins)
 10   continue
      bins=nbins
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindt(npt,time,mag,merr,tbin,nlow,nhigh)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
c     bins data into equal time bins
      implicit none
      integer npt,i,nbins,bin,nmax,j,ntemp,ntemp2,nlow,nhigh
      parameter(nmax=600000)
      real time(npt),mag(npt),merr(npt),tbin,tmin,tmax,ltime,
     .     avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .     var(nmax),ep(nmax),s,p,avgs(nmax),sigcut,xt(nmax),yt(nmax),
     .     zt(nmax),xt2(nmax),yt2(nmax),zt2(nmax)


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

C     reject high and low values
      ntemp2=0
      if((nlow.gt.0).or.(nhigh.gt.0)) then
C     this will probably be really slow.. ugh.
         do 30 i=1,nbins
            ntemp=0
            do 31 j=1,npt
               bin=int(real(nbins)*(time(j)-tmin)/ltime)+1
               if(bin.eq.i) then
                  ntemp=ntemp+1
                  xt(ntemp)=time(j)
                  yt(ntemp)=mag(j)
                  zt(ntemp)=merr(j)
               endif
 31         continue
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
      subroutine bindt2(npt,time,mag,merr,tbin,nlow,nhigh)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
c     bins data into equal time bins
      implicit none
      integer npt,i,nbins,bin,nmax,j,ntemp,ntemp2,nlow,nhigh
      parameter(nmax=600000)
      real time(npt),mag(npt),merr(npt),tbin,tmin,tmax,ltime,
     .     avgm(nmax),avgt(nmax),avge(nmax),abin(nmax),stdev(nmax),
     .     var(nmax),ep(nmax),s,p,avgs(nmax),sigcut,xt(nmax),yt(nmax),
     .     zt(nmax),xt2(nmax),yt2(nmax),zt2(nmax)


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
      nbins=int(ltime/tbin+0.5)+1
      write(6,*) "tbin:",ltime,nbins,tmin,tmax

C     reject high and low values
      ntemp2=0
      if((nlow.gt.0).or.(nhigh.gt.0)) then
C     this will probably be really slow.. ugh.
         do 30 i=1,nbins
            ntemp=0
            do 31 j=1,npt
               bin=int(real(nbins)*(time(j)-tmin)/ltime)+1
c               write(6,*) "bin:",bin,time(j),mag(j),merr(j)
               if(bin.eq.i) then
                  ntemp=ntemp+1
                  xt(ntemp)=time(j)
                  yt(ntemp)=mag(j)
                  zt(ntemp)=merr(j)
               endif
 31         continue
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
c         write(6,*) "avg:",avgs(i),stdev(i)
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rejhilow(npt,x,y,z,nlow,nhigh)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
      implicit none
      integer npt,nlow,nhigh,ntemp,nmax,i
      parameter(nmax=600000)
      real x(npt),y(npt),z(npt),xt(nmax),yt(nmax),zt(nmax)

      call sort(npt,y,x,z)

      ntemp=0

      do 10 i=1+nlow,npt-nhigh
         ntemp=ntemp+1
         xt(ntemp)=x(i)
         yt(ntemp)=y(i)
         zt(ntemp)=z(i)
 10   continue
      
      do 11 i=1,ntemp
         x(i)=xt(i)
         y(i)=yt(i)
         z(i)=zt(i)
 11   continue
      npt=ntemp

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

C      s=0.
C      do 11 i=1,npt
C         s=s+pts(i)
C 11   continue
C      mean=s/npt

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
      subroutine sort(n,x,y,z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,nmax
      parameter(nmax=600000)
      integer p(nmax)
      real x(n),y(n),z(n),tx(nmax),ty(nmax),tz(nmax)

      call rqsort(n,x,p)

      do 10 i=1,n
         tx(i)=x(p(i))
         ty(i)=y(p(i))
         tz(i)=z(p(i))
 10   continue
      
      do 15 i=1,n
         x(i)=tx(i)
         y(i)=ty(i)
         z(i)=tz(i)
 15   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sort2(n,x,y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,nmax
      parameter(nmax=600000)
      integer p(nmax)
      real x(n),y(n),tx(nmax),ty(nmax)

      call rqsort(n,x,p)

      do 10 i=1,n
         tx(i)=x(p(i))
         ty(i)=y(p(i))
 10   continue
      
      do 15 i=1,n
         x(i)=tx(i)
         y(i)=ty(i)
 15   continue

      return
      end


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
