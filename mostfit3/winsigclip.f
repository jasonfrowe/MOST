      program windowsigmaclipping
      implicit none
      integer npt,nmax,nunit,i,iter,iargc
      parameter(nmax=400000)
      real rmavg,time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),
     .   yc(nmax),fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),skystd(nmax),
     .   maglo,maghi,avg,var
      double precision dtime(nmax)
      character*80 filename,cline
      
      iter=1
      
c      write(6,*) "Enter MOST photometry filename"
c      read(5,*) filename
      
      nunit=10 !set unit number for file read

      if(iargc().lt.1)goto 902
      call getarg(1,filename)

      maglo=-1.0
      maghi=-1.0
      
      if(iargc().ge.2)then
        call getarg(2,cline)
        read(cline,*) maglo
      endif
      
      if(iargc().ge.3)then
        call getarg(3,cline)
        read(cline,*) maghi
      endif
      
      if(iargc().ge.4)then
         call getarg(4,cline)
         read(cline,*) iter
      endif
      
      
      open(unit=nunit,file=filename,status='old',err=901)
      
C     read in photometry data  
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
C     we have photometry, can close the file now.
      close(nunit)
      
      if(maglo.eq.maghi)then
        call avevar(mag,npt,avg,var)
        maglo=avg-4.0*sqrt(var)
        maghi=avg+4.0*sqrt(var)
      endif
      write(6,*) "mag-hi-lo:",maglo,maghi
      
      do 10 i=1,iter
         call winsigclip(npt,time,dtime,mag,merr,sky,xc,yc,
     .      fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,
     .      maglo,maghi)
 10   continue
      
      filename="sigclipped.dat"
      rmavg=0.
      call exportdata(npt,dtime,mag,merr,filename,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
      
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 902  write(6,*) "Usage: winsigclip filename maglo maghi <iter>"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine winsigclip(npt,time,dtime,mag,merr,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,
     .     maglo,maghi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,nwidth,i,j,k
      parameter(nmax=400000)
      integer cut(nmax),iord(nmax)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),
     .   xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),
     .   mfield(nmax),aflux(nmax),nap(nmax),itime(nmax),skystd(nmax),
     .   mean,ty(nmax),std,stdev,sigcut,ftotal(nmax),maglo,maghi
      double precision dtime(nmax)

c      maglo=14.6
c      maghi=14.8
c      maglo=16.8
c      maghi=17.4

      sigcut=3.0 !sigma cutting factor
      nwidth=50  !amount of data to sigclip
C     first we sort the data
      call rqsort(npt,time,iord)
      
      do 10 i=1,npt,nwidth    !loop over all data
         k=0  !initialize counter
         mean=0.  !initialize mean
         do 20 j=i,nwidth+i-1 !loop over subset
            if((j.lt.npt).and.(mag(iord(j)).gt.maglo).and.
     .           (mag(iord(j)).lt.maghi)) then
               k=k+1
               mean=mean+mag(iord(j))
               ty(k)=mag(iord(j))
            endif
 20      continue
         mean=mean/real(k)
 21      std=stdev(k,ty,mean)
         do 30 j=i,nwidth+i-1
            if(j.lt.npt) then
               if(abs(mag(iord(j))-mean).lt.sigcut*std) then
                  cut(iord(j))=0 !data not clipped
               else
                  cut(iord(j))=1 !data clipped
               endif
c               if(fy(iord(j)).gt.1.25) cut(iord(j))=1
c               if(sky(iord(j)).lt.400.0) cut(iord(j))=1
c               if(sky(iord(j)).gt.650.0) cut(iord(j))=1
            endif
 30      continue
 10   continue
      
      j=0 !initialize counter
      do 40 i=1,npt
         if(cut(i).eq.0)then
            j=j+1
            time(j)=time(i)
            dtime(j)=dtime(i)
            mag(j)=mag(i)
            merr(j)=merr(i)
            sky(j)=sky(i)
            xc(j)=xc(i)
            yc(j)=yc(i)
            fx(j)=fx(i)
            fy(j)=fy(i)
            fxy(j)=fxy(i)
            tboard(j)=tboard(i)
            mfield(j)=mfield(i)
            ftotal(j)=ftotal(i)
            aflux(j)=aflux(i)
            nap(j)=nap(i)
            itime(j)=itime(i)
            skystd(j)=skystd(i)
         endif
 40   continue   
      npt=j
     
 999  return
      end
      
      
     
