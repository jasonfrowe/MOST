      program boothisto
      implicit none
      integer nmax,npt,col,nbin,nfreq,i,dumi,j,k
      parameter(nmax=500000)
      real dd(nmax),xb1,xb2,yb1,yb2,datamin,datamax,bdatax(nmax),
     .   bdatay(nmax),bmax,fc,ave,var,pi,dif,ans(2,nmax),frsol(nmax),
     .	 errs(6)
      character*80 filename,freqname,cline
      pi=acos(-1.0)
      
      nbin=30 !default binning for histograms
      if(iargc().lt.2) goto 903 !need two arguements for this program.
      call getarg(1,freqname)
      open(unit=11,file=freqname,status='old',err=901)
      call getarg(2,filename)
      open(unit=10,file=filename,status='old',err=901)
      if(iargc().ge.3)then
      	call getarg(3,cline)
      	read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 902
      
C     First we read in the Frequency Analysis data.
C     This tells us the number of frequencies      
      call readfreqs(nfreq,nmax,frsol)
      close(11) !close the file, we don't need it anymore.
      
      call pgopen('?') !open PGPlot device
      call pgask(.false.) !don't ask for new page.. just do it.
      call pgpage() !start a multi-page document
C     Now we loop over all Frequencies.
      do 10 i=1,3*nfreq !there are 3 variables for each frequency.
         col=i
         npt=nmax
         call getdata(col,npt,dd) !read column of bootstrapped values

C     These lines will sigclip at 4 sigma if required.  (probably not)
c         call avevar(dd,npt,ave,var)
c         k=0
c         do 12 j=1,npt
c            if(abs(dd(j)-ave).lt.4.0*sqrt(var))then
c               k=k+1
c               dd(k)=dd(j)
c            endif
c 12      continue
c         npt=k
C     End of sigmaclipping

         call avevar(dd,npt,ave,var)
C     If we have a phase parameter, then we need to check for cyclic 
C     numbers around Pi
         if(mod(i,3).eq.0)then !see if this is the 3rd variable
            dif=ave-pi !center on Pi
            do 11 j=1,npt
C     now we shift all data around the midpoint
               if(dd(j)-dif.gt.2.0*Pi)dd(j)=dd(j)-2.0*Pi
               if(dd(j)-dif.lt.0.0)dd(j)=dd(j)+2.0*Pi 
 11         continue
            call avevar(dd,npt,ave,var)
         endif
         ans(1,i)=ave
         ans(2,i)=sqrt(var)
      
         call findrange(npt,dd,datamin,datamax)
c         write(6,*) "datamin,datamax",datamin,datamax
         call bindata(nbin,npt,dd,bdatax,bdatay,datamin,datamax,bmax)
      
         call pgpage() !fresh plotting surface
         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
         call pgvport(xb1,xb2,yb1,yb2)
         fc=2.*(datamax-datamin)/real(nbin) !center histograms
         call pgwindow(datamin+fc,datamax+fc,0.,bmax+0.1*bmax) !set size
         call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
         call pglabel("Value","Frequency","") !add axis labels
         call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
          
         call errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol(i),errs)

		 if(mod(i+1,3).eq.0)then
		 	write(6,500) 1000.0*frsol(i),1000.0*ans(1,i),
     .			1000.0*ans(2,i),(1000.0*errs(j),j=1,6)
     	 else
         	write(6,500) frsol(i),ans(1,i),ans(2,i),(errs(j),j=1,6)
		 endif
		 if(mod(i,3).eq.0)write(6,*)
 500     format(9(F11.7,1X))
C 500     format(9(F15.12,1X))
         
C        this line is heavy on disk-io but easy on memory.
         rewind(10) ! rewind the bootstrap file so we can read it again.
         read(10,*) dumi !not needed, but no harm. Just skips first line
 10   continue
      call pgclos() !close plotting device
      
c      do 12 i=0,nfreq-1
c         write(6,500) (ans(1,j),j=3*i+1,3*i+3),(ans(2,j),j=3*i+1,3*i+3)
c 12   continue
c 500  format(200(E14.7,1X))
      
      
      close(10)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 902  write(6,*) "Error: Nbin must be at least 2"
      goto 999
 903  write(6,*) "Usage: boothist <freqname> <bootstats> [nbin]"
      write(6,*) "<freqname> : Frequency output list from MOSTPer"
      write(6,*) "<bootstats>: Bootstrap statistics from COSBoot"
      write(6,*) "[nbin] (optional) : binned parameter for histrograms, 
     .default=30"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,nfreq,i,j,mdpnt,k
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .	intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
      	if((bdatax(i).lt.frsol).and.(bdatax(i+1).gt.frsol))then
        	mdpnt=i
      	endif
 5    continue
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initalize variable to zero.
        errs(i)=0.0
 14	  continue
      do 13 k=1,bmax-1
      	shgt=real(bmax-k)
C     	first we move in the forward direction      
      	loop=.true. !loop until loop is false
      	i=mdpnt-1
      	do 10 while(loop)
      		if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
      			loop=.false.
c      			inth=bdatax(i)
      			inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .				*(bdatax(i+1)-bdatax(i))+
     .				bdatax(i)
      		endif
        	if(i.ge.nbin-1) then
        		loop=.false.
        		inth=bdatax(nbin)
    		endif
    		i=i+1
c    		write(6,*) "i:",i
 10   	enddo
        
C     	now we move in the reverse direction
      	loop=.true. !loop until loop is false
		i=mdpnt+1
      	do 11 while(loop)
      		if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
      			loop=.false.
c      			intl=bdatax(i)
				intl=bdatax(i)-
     .				(bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .				*(bdatax(i)-bdatax(i-1))
      		endif
        	if(i.le.2) then
        		loop=.false.
        		intl=bdatax(1)
    		endif
    		i=i-1
 11   	enddo

 	  	j=0
      	do 12 i=1,npt
      		if((dd(i).gt.intl).and.(dd(i).lt.inth))then
      			j=j+1    
      		endif
 12	  	continue
 
 		per=real(j)/real(npt)
 		if((per.ge.0.683).and.(perold.lt.0.683))then
      		e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
      		errs(1)=e1-frsol
      		errs(2)=e2-frsol
c			write(6,*) "1 sigma:",frsol,errs(1),errs(2)
 		endif
 		if((per.ge.0.954).and.(perold.lt.0.954))then
      		e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
      		errs(3)=e1-frsol
      		errs(4)=e2-frsol
c			write(6,*) "2 sigma:",frsol,errs(3),errs(4)
 		endif 
 		if((per.ge.0.9973).and.(perold.lt.0.9973))then
      		e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
      		errs(5)=e1-frsol
      		errs(6)=e2-frsol
c			write(6,*) "3 sigma:",frsol,errs(5),errs(6)
 		endif     	
      	perold=per
      	intlold=intl
      	inthold=inth
c      	write(6,*) per,intl,inth,shgt
 13	  continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nfreq,nmax,frsol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfreq,nmax,nvar,i,j
      real frsol(nmax)
      double precision ztime
      
      read(11,*) ztime
      read(11,*) nvar,nfreq
      
      do 10 i=1,nfreq
        j=3*(i-1)+1 !read all variables in order.
      	read(11,*) frsol(j),frsol(j+1),frsol(j+2)
 10   continue
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
      ave=0.0
      do 11 j=1,n
         ave=ave+data(j)
 11   continue
      ave=ave/n
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
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.lt.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i)*binsize
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(npt,pdata,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real pdata(npt),datamin,datamax
      
C     initialize datamin and datamax
      datamin=pdata(1)
      datamax=pdata(1)
      
      do 10 i=2,npt
         datamin=min(pdata(i),datamin)
         datamax=max(pdata(i),datamax)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowsetup(xb1,xb2,yb1,yb2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      setup plotting area
      real xb1,xb2,yb1,yb2,xs1,xs2,
     .   ys1,ys2,pratio
      
C     set text size      
      call pgsch(1.0)
C     inquire window size
      call pgqvsz(2,xs1,xs2,ys1,ys2)
      pratio=(xs2-xs1)/(ys2-ys1)
C     if we have a landscape surface
      if(xs2-xs1.gt.ys2-ys1)then
C        start 10% from the edge
         yb1=(ys2-ys1)*0.10+ys1
C        use 60% of surface available
         yb2=(ys2-ys1)*0.70+yb1
C        now normalize to 1
         yb1=yb1/(ys2-ys1)         
         yb2=yb2/(ys2-ys1)
C        now use same size in other dimension      
         xb1=yb1
         xb2=yb1+(yb2-yb1)/pratio
      else
C        start 10% from the edge
         xb1=(xs2-xs1)*0.10+xs1
C        use 60% of surface available
         xb2=(xs2-xs1)*0.60+xb1
C        now normalize to 1
         xb1=xb1/(xs2-xs1)         
         xb2=xb2/(xs2-xs1)
C        now use same size in other dimension      
         yb1=xb1
         yb2=xb1+(xb2-xb1)*pratio
      endif
      
      return
      end 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(col,npt,dd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer col,npt,i,j
      real dd(npt),dumr
      
      i=1
 10   read(10,*,end=11,err=901) (dumr,j=1,col-1),dd(i)
c         write(6,*) dd(i)
         i=i+1
         goto 10
 11   continue
      
      npt=i-1
      
 500  format(200(E14.7,1X))
      goto 999
 901  write(6,*) "Error on line: ",i
      goto 999
      
 999  return
      end

