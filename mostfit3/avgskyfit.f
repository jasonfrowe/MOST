      program avgskyfit
      implicit none
      integer nmax,nunit,npt,nbin,i,nptb,sp
      parameter (nmax=400000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   aflux(nmax),nap(nmax),itime(nmax),ftotal(nmax),skybin(nmax),
     .   magbin(nmax),rmavg,maglow,maghi,btime(nmax),bmag(nmax),
     .   bmerr(nmax),bsky(nmax),stdsky(nmax),bstdsky(nmax)
      double precision dtime(nmax)
      character*80 filename,outname
      
      sp=2
      maglow=-0.4
      maghi=0.1
      
      nunit=10  !set unit number for file.
      write(6,*) "Enter Photometry File "
      read(5,*) filename !read in filename
    
C     open up filename for reading.  File must exist or exit with
C     an error code  
      open(unit=nunit,file=filename,status='old',err=901)
     
      write(6,*) "Enter Output Photometry name "
      read(5,*) outname
     
C     call subroutine to read in data. 
      call readdataskyfit(nunit,npt,time,dtime,mag,merr,
     .     sky,xc,yc,fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,
     .     stdsky,maglow,maghi)

      call copybin(npt,time,mag,merr,sky,stdsky,mfield,nptb,btime,bmag,
     .   bmerr,bsky,bstdsky,sp)

      call makeskybin(npt,mag,merr,sky,stdsky,nptb,btime,bmag,bmerr,
     .   bsky,bstdsky,nbin,magbin)
      
      do 10 i=1,npt
         mag(i)=mag(i)-magbin(i)
 10   continue
      
      rmavg=0.0
      call exportdata(npt,dtime,mag,merr,outname,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,stdsky,rmavg)
     
      close(nunit)
      goto 999
 901  write(6,*) "Error: Cannot open ",filename
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine makeskybin(npt,mag,merr,sky,stdsky,nptb,btime,bmag,
     .   bmerr,bsky,bstdsky,nbin,magbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,nmax,j,k,nptb
      parameter(nmax=400000)
      real mag(npt),merr(npt),sky(npt),magbin(npt),cmean,
     .   binwidth,tempmag(nmax),tempsky(nmax),tempmerr(nmax),sig,
     .   btime(nptb),bmag(nptb),bmerr(nptb),bsky(nptb),bstdsky(nptb),
     .   stdsky(nmax)
     
      sig=3.0
      
      do 10 i=1,npt
C         binwidth=log10(10.0*sky(i))*10.0+10.0
co         binwidth=log10(sky(i))         
         binwidth=log10(stdsky(i)*stdsky(i))
         k=0
         do 11 j=1,nptb
co            if(abs(bsky(j)-sky(i)).lt.binwidth)then
            if(abs(bstdsky(j)-stdsky(i)).lt.binwidth)then
               k=k+1
c               tempsky(k)=bsky(j)
               tempsky(k)=bstdsky(j)
               tempmag(k)=bmag(j)
               tempmerr(k)=bmerr(j)
            endif
 11      continue
c         write(6,*) i,k
         call sigclip(k,tempsky,tempmag,tempmerr,sig)
c         write(6,*) "hello",k
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
c      mfield(i)=25000.0
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
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine copybin(npt,time,mag,merr,sky,stdsky,mfield,nptb,btime,
     .   bmag,bmerr,bsky,bstdsky,sp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nptb,i,sp
      real time(npt),mag(npt),merr(npt),btime(npt),bmag(npt),bmerr(npt),
     .   mfield(npt),sig,x1,x2,sky(npt),bsky(npt),stdsky(npt),
     .   bstdsky(npt)
      
            nptb=0
      do 5 i=1,npt
c         if(mfield(i).gt.20000.0) then
            if(sp.eq.2) then
               if(sky(i).lt.1.0) goto 60
               x1=2042.83
               x2=2043.13
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1-3.52/2.0
               x2=x2-3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52
               x2=x2+3.52
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2046.35
               x2=2046.65
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2049.86
               x2=2050.16
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2053.39
               x2=2053.69
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2056.91
               x2=2057.21
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2060.44
               x2=2060.74
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2063.98
               x2=2064.28
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2067.49
               x2=2067.79
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2071.01
               x2=2071.31
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2074.54
               x2=2074.84
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2078.08
               x2=2078.38
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=2081.60
               x2=2081.90
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
            elseif(sp.eq.1) then
               x1=1690.4
               x2=1690.6
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=1693.8
               x2=1694.2
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=1697.4
               x2=1697.7
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=1700.9
               x2=1701.2
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
               x1=x1+3.52/2.0
               x2=x2+3.52/2.0
               if((time(i).gt.x1).and.(time(i).lt.x2)) goto 60
            endif  
            nptb=nptb+1
            btime(nptb)=time(i)
            bmag(nptb)=mag(i)
            bmerr(nptb)=merr(i)
            bsky(nptb)=sky(i)
            bstdsky(nptb)=stdsky(i)
 60      continue
c         endif
 5    continue
c      sig=3.0
c      do 6 i=1,2
c         call sigclip(nptb,btime,bmag,bmerr,sig)
c 6    continue
      write(6,*) "nptb1",nptb
      
      return
      end