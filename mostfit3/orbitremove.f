c23456789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7
      program orbitremove
      implicit none
      integer nmax,nunit,npt,nphase,nfold,nptb,i,sp
      parameter(nmax=500000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),period,
     .   bphase(nmax),mindate,btime(nmax),bmag(nmax),bmerr(nmax),tbin,
     .   rmavg,sig,stdsky(nmax)
      double precision dtime(nmax)
      character*80 filename,ans
      
C     Get the filename from the user
      write(6,*) "Enter MOST photometry data file"
      read(5,500) filename
 500  format(A80)
      
C     Assign a unit number and then attempt to open the file
C     for access, if we encounter an error we jump to the end of the 
C     program were error messages are parsed.
      nunit=10
      open(unit=nunit,file=filename,status='old',err=901)

C     Standard subroutine for reading in data.
      write(6,*) "Reading in data"
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,fy,fxy,
     .   tboard,mfield,ftotal,aflux,nap,itime,stdsky)
     
 11   write(6,*) "Date type: 0,1,2 [0]"
      read(5,500) ans
      if(ans.eq." ") then
         sp=0
      else
         read(ans,*,err=11) sp
      endif
      if((sp.lt.0).or.(sp.gt.2)) then
         write(6,*) "That makes no sense"
         goto 11
      endif
           
C     close input file
      close(10)

C     make a copy of data for binning
      call copybin(npt,time,mag,merr,sky,mfield,nptb,btime,bmag,bmerr,
     .   sp)

C     bin the data, tbin is in minutes
 12   write(6,*) "Enter bin interval (in minutes) [0.0]"
      read(5,500) ans
      if(ans.eq." ") then
         tbin=0.0
      else
         read(ans,*,err=12) tbin
      endif
      if(tbin.lt.0.0)then
         write(6,*) "That makes no sense"
         goto 12
      endif
      
      if(tbin.gt.0.0) then
         call bindt(nptb,btime,bmag,bmerr,tbin,0,0)
         write(6,*) "nptb",nptb
      else
         do 13 i=1,npt
            btime(i)=time(i)
            bmag(i)=mag(i)
            bmerr(i)=merr(i)
 13      continue
         nptb=npt   
      endif

 14   write(6,*) "Enter Period to fold data with [MOSTorb]"
      read(5,500) ans
      if(ans.eq." ") then
         period=0.0704256
      else
         read(ans,*,err=14) period
      endif
      if(period.lt.0.) then
         write(6,*) "Answer makes no sense"
         goto 14
      endif

C     Okay.. lets start by phasing the data to the orbital period.
C     using double precision to preserve time accuracy
      call rphase(nptb,btime,bphase,period,nphase,mindate)
     
C     Report number of phases found
      write(6,501) "Number of phases found:",nphase     
 501  format(A23,1X,I4)
     
C     Get number of phases to fold with
 10   write(6,502) "How manys to fold together? [",nphase,"]"
 502  format(A29,I5,A1)
      read(5,500) ans
      if(ans.eq." ") then
         nfold=nphase
      else
         read(ans,*,err=10) nfold
      endif
c      if((nfold.le.0).or.(nfold.gt.nphase))then
c         write(6,*) "Answer makes no sense.. try again"
c         goto 10
c      endif

C     alright... lets do it.
      call foldnremove(npt,time,mag,nptb,bphase,bmag,bmerr,
     .   nfold,mindate,period)

      rmavg=0.
      filename="newdata.dat"
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,stdsky,rmavg)
           
C     Error calls are here
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine copybin(npt,time,mag,merr,sky,mfield,nptb,btime,bmag,
     .   bmerr,sp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nptb,i,sp
      real time(npt),mag(npt),merr(npt),btime(npt),bmag(npt),bmerr(npt),
     .   mfield(npt),sig,x1,x2,sky(npt)
      
            nptb=0
      do 5 i=1,npt
         if(mfield(i).gt.20000.0) then
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
            
            
 60      continue
         endif
 5    continue
      sig=3.0
      do 6 i=1,0
         call sigclip(nptb,btime,bmag,bmerr,sig)
 6    continue
      write(6,*) "nptb1",nptb
      
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine foldnremove(npt,time,mag,nptb,bphase,bmag,
     .   bmerr,nfold,mindate,period)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none 
      integer npt,nfold,i,j,nptb,nmax,nb,n
      parameter(nmax=500000)
      integer nd(nmax)
      real time(npt),mag(npt),bphase(nptb),bmag(nptb),mindate,yp1,yp2,
     .   y2(nmax),bx(nmax),by(nmax),phase,pshift,y,period,perc,percold,
     .   ymin,ymax,avgx(nmax),avgy(nmax),avgz(nmax),bz(nmax),sig,
     .   bmerr(nmax),pbin,avgbin
      character*80 tmpc,ans
      
      pbin=0.005
 5    write(6,*) "Enter phase size for binning [0.005]"
      read(5,501) ans
 501  format(A80)
      if(ans.eq." ") then
         pbin=0.005
      else
         read(ans,*,err=5) pbin
      endif
      if((pbin.le.0.0).or.(pbin.ge.1.0))then
         write(6,*) "That answer makes no sense"
         goto 5
      endif
             
      
      percold=0.
      
c      call pgopen('?')
      
      open(unit=12,file="shaperm.dat")
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
         avgbin=0.
         do 14 j=1,nb
            avgbin=avgbin+by(j)
            if(abs(0.5-bx(j)).lt.pbin) then
               y=y+by(j)
               n=n+1
               avgx(n)=bx(j)
               avgy(n)=by(j)
               avgz(n)=bz(j)
            endif
 14      continue
c         avgbin=avgbin/real(nb)
          call rqsort(nb,by,nd)
          avgbin=by(nd(nb/2))
         
         sig=3.0
         call sigclip(n,avgx,avgy,avgz,sig)
         
c         y=0
c         do 15 j=1,n
c            y=y+avgy(j)
c 15      continue
c         y=y/real(n)
C     Find medians
      call rqsort(n,avgy,nd)
      y=avgy(nd(n/2))

c         call pgsci(3)
c         call pgpt1(0.5,y,10)
c         call pgsci(1)
         mag(i)=mag(i)-y+avgbin
         write(12,502) time(i),mag(i),y
 502     format(F13.8,1X,2(F9.6,1X))
 10   continue
      close(12)
 
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
      mindate=time(1)
      do 5 i=2,npt
         mindate=min(mindate,time(i))
 5    continue
      
C     find phase of date (not integer chopping)

      nphase=int((time(1)-mindate)/period)
      do 10 i=2,npt
         phase(i)=(time(i)-mindate)/period
C        nphase is the larger integer phase found
         nphase=max(nphase,int(phase(i)))
 10   continue
 
      return
      end
      