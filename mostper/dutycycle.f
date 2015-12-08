      program dutycycle
      implicit none
      integer nmax,i,nunit,npt,j
      parameter(nmax=500000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),sig,
     .   itime(nmax),skystd(nmax),dc
      double precision dtime(nmax),deltat(nmax),dtot
      character*80 filename
      
      write(6,*) "Enter MOST data file name"
      read(5,*) filename
      
      nunit=10
      open(unit=nunit,file=filename,status="old",err=901)
      
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
      
      sig=3.0
      do 11 i=1,0
         call sigclip(npt,dtime,time,mag,merr,sig)
 11   continue      
      
      j=0
      dtot=0.0d0
      do 10 i=1,npt-1
         deltat(i)=dtime(i+1)-dtime(i)
         deltat(i)=deltat(i)*24.0*60.0*60.0
c         write(6,*) deltat(i)
c         read(5,*)
         if(deltat(i).lt.85.0) then
            j=j+1
            dtot=dtot+deltat(i)
         endif
 10   continue
      dc=dtot/((dtime(npt)-dtime(1))*24.0*60.0*60.0)
      write(6,*) "Duty Cycle: ",dc
         
      
      close(10)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigclip(npt,dtime,time,mag,merr,sig)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=500000,nitmax=50)
      real time(npt),mag(npt),merr(npt),sig,std,stdev,tmp1(nmax),
     .     tmp2(nmax),tmp3(nmax),mean,dum,omean
      double precision dtime(npt),tmp4(nmax)

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
            tmp4(ntmp)=dtime(i)
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
         dtime(i)=tmp4(i)
 20   continue
      npt=ntmp

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