      program frequencyremover
      implicit none
      integer addsub,nunit,npt,nmax,nfreq,nunit2,iargc,i
      parameter(nmax=500000)
      real time(nmax),mag(nmax),merr(nmax),sky(nmax),xc(nmax),yc(nmax),
     .   fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),mfield(nmax),
     .   ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),skystd(nmax),
     .   freq(nmax),amp(nmax),ph(nmax),ztime,rmavg
      double precision dtime(nmax)
      character*80 filename,cline,outname,freqname
      
      if(iargc().lt.4)goto 900
      call getarg(1,filename)
      call getarg(2,cline)
      read(cline,*) addsub
      if((addsub.lt.0).or.(addsub.gt.1))goto 900
      call getarg(3,freqname)
      call getarg(4,outname)
      
c 10   write(6,*) "Add[0] or remove[1] frequencies?"
c      read(5,*,err=10) addsub
c      if((addsub.lt.0).or.(addsub.gt.1))goto 10 
      
c      write(6,*) "Enter MOST photometry filename"
c      read(5,*) filename
      
      nunit=10 !set unit number for file read
      
      open(unit=nunit,file=filename,status='old',err=901)
      
c      write(6,*) "Enter file containing frequency information"
c      read(5,*) freqname
      
C     open file with frequencies, amplitudes and phases
C     from a cosine series
      nunit2=12
      open(unit=nunit2,file=freqname,status='old',err=901)
    
C     read in photometry data  
      call readdata(nunit,npt,time,dtime,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
C     we have photometry, can close the file now.
      close(nunit)
      
C     These lines convert to mmag
      do 11 i=1,npt
        mag(i)=mag(i)*1000.0
        merr(i)=merr(i)*1000.0
 11   continue   
      
C     read in frequency information.
      nfreq=nmax
      call readfreqs(nunit2,nfreq,ztime,freq,amp,ph)
C     have frequencies, can close file now.
      close(nunit2)
      
      call freqaddsub(addsub,npt,nfreq,ztime,time,mag,freq,amp,ph)

C     These lines convert to mag
      do 12 i=1,npt
        mag(i)=mag(i)/1000.0
        merr(i)=merr(i)/1000.0
 12   continue         
      
      
      rmavg=0.
      call exportdata(npt,dtime,mag,merr,outname,sky,xc,yc,
     .     fx,fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd,rmavg)
      
      goto 999
 900  write(6,*) "Usage: freqremove <fname> <addsub> <freqdat> <outdat>"
		write(6,*) "<fname>: MOST Photometry datafile"
		write(6,*) "<addsub>: 0 - add frequencies"
		write(6,*) "          1 - remove frequencies"
		write(6,*) "<freqdat>: Frequency output from MOSTPer"
		write(6,*) "<outdat>: Output name for photometry"
		goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine freqaddsub(addsub,npt,nfreq,ztime,time,mag,freq,amp,ph)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer addsub,npt,nfreq,i,j
      real time(npt),mag(npt),freq(nfreq),amp(nfreq),ph(nfreq),Pi,tPi,
     .   ztime,ttime,tmag
      Pi=3.141592654
      tPi=2.0*Pi
      
      
      do 10 i=1,npt
         ttime=time(i)-ztime
         tmag=0.
         do 11 j=1,nfreq
            tmag=tmag+amp(j)*cos(tPi*freq(j)*ttime+ph(j))
 11      continue
         if(addsub.eq.0) then
            mag(i)=mag(i)+tmag
         elseif(addsub.eq.1) then
            mag(i)=mag(i)-tmag
         endif
 10   continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nunit,nfreq,ztime,freq,amp,ph)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nfreq,i,dumi
      real freq(nfreq),amp(nfreq),ph(nfreq),ztime
      
      
      read(nunit,*) ztime
      read(nunit,*) dumi,dumi
      i=1
 10   read(nunit,*,end=11) freq(i),amp(i),ph(i)
         i=i+1
         goto 10
 11   continue
      nfreq=i-1
      
      return
      end
      
      