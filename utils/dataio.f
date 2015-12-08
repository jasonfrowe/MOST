C Contents
C
C readdata
C ntest
C ovrwrt

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in data
C     npt - number of data points (returned)
C     id - id of star (returned)
C     xcoo,ycoo - co-ordinates of star (returned)
C     time,mag,merr - the data and associated error.
      implicit none
      integer nmax,nunit
      parameter(nmax=600000)
      integer npt,i,ntest,rn,filet
      real time(nmax),mag(nmax),merr(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),skystd(nmax),g,Ist,r,dumr
      double precision dtime(nmax)
      character*10 dumc

      i=1

      filet=0

      if(filet.eq.0)then
         goto 10
   9     continue
            write(6,*) "error reading",i
  10        read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i),sky(i),
     .         xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .         ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
c     .         ftotal(i),dumc,nap(i),itime(i),skystd(i)
c            read(dumc,*) aflux(i)
            time(i)=real(dtime(i))
            merr(i)=abs(merr(i))
            if(merr(i).le.0.0) goto 10
c        if(fy(i).lt.0.5) goto 10
c        if(fx(i).lt.0.5) goto 10
c        if((mag(i).lt.12.15).or.(mag(i).gt.12.45)) goto 10
C     Check for NaNs
            rn=ntest(mag(i))
            if(rn.eq.1) goto 10
            rn=ntest(sky(i))
            if(rn.eq.1) goto 10
            i=i+1
            goto 10
 20      continue
      elseif(filet.eq.1)then
 29      continue
 30      read(10,*,err=29,end=40) dtime(i),mag(i),merr(i),sky(i),dumr,
     .      dumr,dumr,dumr,itime(i)
           xc(i)=0.
           yc(i)=0.
           fx(i)=0.
           fy(i)=0.
           fxy(i)=0.
           ftotal(i)=0.
           tboard(i)=0.
           mfield(i)=25000.0
           aflux(i)=0.
           nap(i)=0.
           skystd(i)=0.
           time(i)=dtime(i)
c           if(sky(i).gt.2000.0) goto 30
           if(mag(i).gt.80.0) goto 30
           merr(i)=abs(merr(i))
           if(merr(i).le.0.0) goto 30
           i=i+1
         goto 30
 40      continue
      endif
      
C     old data format.
c 505  format(F10.5,1X,3(F7.4,1X),F8.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
c     .     1X,F7.3,1X,F9.2)
C             13.8/11.6

 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
      npt=i-1
      write(6,*) "Points read: ",npt

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real mag(npt),merr(npt),sky(npt),xc(npt),
     .     yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),mfield(npt),
     .     rmavg,ftotal(npt),aflux(npt),nap(npt),itime(npt),skystd(npt)
      double precision dtime(npt)
      character*80 filename

      open(unit=11,file=filename)

      do 10 i=1,npt
         if(mag(i)+rmavg.lt.80.0) then
            write(11,510) dtime(i),mag(i)+rmavg,merr(i),sky(i),
     .           xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .         ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
         endif
 10   continue

 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
      close(11)

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
