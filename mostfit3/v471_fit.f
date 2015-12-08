      program v471tau
      implicit none
      integer nmax,nunit,npt,nbin,i
      parameter (nmax=400000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   aflux(nmax),nap(nmax),itime(nmax),ftotal(nmax),skybin(nmax),
     .   magbin(nmax),rmavg,skystd(nmax)
      double precision dtime(nmax)
      character*80 filename,outname
      
      nunit=10  !set unit number for file.
      write(6,*) "Enter V471 Tau Photometry File "
      read(5,*) filename !read in filename
    
C     open up filename for reading.  File must exist or exit with
C     an error code  
      open(unit=nunit,file=filename,status='old',err=901)
     
      write(6,*) "Enter Output Photometry name "
      read(5,*) outname
     
C     call subroutine to read in data. 
      call readdatav471(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,fy,
     .   fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)

      call makeskybin(npt,mag,merr,sky,nbin,skybin,magbin)
      
      do 10 i=1,npt
         mag(i)=mag(i)-magbin(i)
 10   continue
      
      rmavg=0.
      call exportdata(npt,dtime,mag,merr,outname,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
     
      close(nunit)
      goto 999
 901  write(6,*) "Error: Cannot open ",filename
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine makeskybin(npt,mag,merr,sky,nbin,skybin,magbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,nmax,j,k
      parameter(nmax=400000)
      real mag(npt),merr(npt),sky(npt),skybin(npt),magbin(npt),cmean,
     .   binwidth,tempmag(nmax),tempsky(nmax),tempmerr(nmax),sig
     
      sig=3.0
      
      do 10 i=1,npt
         binwidth=log10(10.0*sky(i))*10.0+10.0
         k=0
         do 11 j=1,npt
            if(abs(sky(j)-sky(i)).lt.binwidth)then
               k=k+1
               tempsky(k)=sky(j)
               tempmag(k)=mag(j)
               tempmerr(k)=merr(j)
            endif
 11      continue
         call sigclip(k,tempsky,tempmag,tempmerr,sig)
         skybin(i)=cmean(k,tempsky)
         magbin(i)=cmean(k,tempmag)
c         write(6,*) magbin(i),skybin(i),binwidth
 10   continue      
     
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function cmean(n,f)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer n,i
      real f(n)
      
      cmean=0.
      do 10 i=1,n
         cmean=cmean+f(i)
 10   continue
      cmean=cmean/real(n)
      
      return
      end
      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdatav471(nunit,npt,time,dtime,mag,merr,sky,xc,yc,
     .   fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
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
     .     itime(nmax),skystd(nmax)
      double precision dtime(nmax)

      i=1

   9  continue
  10  read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i),sky(i),
     .     xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),ftotal(i),
     .     aflux(i),nap(i),itime(i),skystd(i)
c      mfield(i)=25000.0
      time(i)=real(dtime(i))
         if(sky(i).gt.14200.0) goto 10
        merr(i)=abs(merr(i))
c  10  read(10,505,err=9,end=20) time(i),amag(i),mag(i),merr(i),sky(i),
c     .     xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i)
        if(merr(i).le.0.0) goto 10
c        if(fy(i).lt.0.5) goto 10
c        if(fx(i).lt.0.5) goto 10
C     1 14.05 14.70
C       13.60 14.20
c        if((mag(i).lt.13.80).or.(mag(i).gt.15.60)) goto 10
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