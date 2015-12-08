      program binplotting
      implicit none
      integer nmax,i,npt,now(3),seed,j,ns,niter,k
      parameter(nmax=400 000,ns=100,niter=100)
      real time(nmax),phase(nmax),mag(nmax),merr(nmax),Eerr,mean,std,
     .   stdev,dumr,ran2,gasdev,rannums(nmax),bprob,mm,mmin,mmax,
     .   rannums2(nmax),plotx(ns),ploty(ns)
      character*80 filename
      
C     Random number creation (so we get a different seed each time)
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)      
      
      write(6,*) "Enter filename :"
      read(5,*) filename
      
      open(unit=10,file=filename,status='old',err=901)
      
      i=1
 10   read(10,*,end=11) time(i),phase(i),mag(i),merr(i)
         i=i+1
         goto 10
 11   continue
      npt=i-1
      close(10)
      
      call pgopen('?')
      
      Eerr=0.
      do 12 i=1,npt
         Eerr=Eerr+merr(i)*merr(i)
 12   continue
      Eerr=Sqrt(Eerr)/Sqrt(real(npt))
      
      std=stdev(npt,mag,mean)
      call plotbins(npt,mag,mean,std,Eerr)
      write(6,*) "Mean :",mean, " Std:",std
    
c     std=stdev(npt,rannums,mean)
c     call plotbins(npt,rannums,mean,std,Eerr)    

      mmin=-5.0e-4
      mmax= 5.0e-4      
      do 16 i=1,ns
         plotx(i)=(mmax-mmin)*real(i)/real(ns)+mmin
         ploty(i)=0.0
 16   continue
 
      do 17 k=1,niter
         do 13 i=1,npt
            rannums(i)=merr(i)*gasdev(seed)
 13      continue  
         do 14 j=1,ns
      
            mm=plotx(j)
         
            do 15 i=1,npt
               rannums2(i)=rannums(i)+mm
 15         continue
         
            call baytest(npt,npt,rannums2,mag,bprob)
            ploty(j)=ploty(j)+bprob
         
 14      continue
 17   continue
 
      do 18 i=1,ns
         ploty(i)=ploty(i)/real(niter)
         write(6,*) plotx(i),ploty(i)
 18   continue
 
      
      
      call pgclos()
      
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine baytest(nrand1,nrand2,rans1,rans2,bprob)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrand1,nrand2,nintg,i
      parameter(nintg=100)
      real rans1(nrand1),rans2(nrand2),bprob,mumax,mumin,mintg,pintg,
     .   xx,yy(nintg),rnintg,buxuy,integrate,plotx(nintg),ploty(nintg),
     .   xb1,xb2,yb1,yb2
      
c      call pgpage()
c      call windowsetup(xb1,xb2,yb1,yb2) !setup a square plot
c      call pgvport(xb1,xb2,yb1,yb2) !make a square plotting surface
c      call pgwindow(-5.0,5.0,0.0,2.0e-4)!set window scale
c      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)!boxes and ticks
c      call pglabel("Mean","Probability","") !label axes
      
      mumin=-20.0 !range of mu to scan
      mumax= 20.0
      
      rnintg=nintg !sort integer into real number
      do 10 i=1,nintg
         xx=mumin+real(i)/rnintg*(-mumin)
         yy(i)=buxuy(nrand1,nrand2,rans1,rans2,xx)!sort probs into array
         plotx(i)=xx
         ploty(i)=yy(i)
 10   continue
c      call pgline(nintg,plotx,ploty)
 
      xx=-mumin/rnintg
      mintg=integrate(nintg,xx,yy) !use Simpson rule to get integral
      
      do 20 i=1,nintg
         xx=real(i-1)/rnintg*mumax
         yy(i)=buxuy(nrand1,nrand2,rans1,rans2,xx)!sort probs into array
         plotx(i)=xx
         ploty(i)=yy(i)
 20   continue
c      call pgline(nintg,plotx,ploty)
 
      xx=mumax/rnintg
      pintg=integrate(nintg,xx,yy) !use Simpson rule to get integral

      bprob=min(mintg/pintg,pintg/mintg)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function buxuy(nptx,npty,datx,daty,muxmy)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nptx,npty,i
      real datx(nptx),daty(npty),muxmy,Sx,Sy,temp(3),avex,avey,v,tp,
     .   rnrx,rnry,s,rgammln,Pi
      Pi=3.141592654
     
      rnrx=Real(nptx) !sort integer into real variable
      rnry=Real(npty)
      
      call avevar(datx,nptx,avex,Sx)  !get average and variance
      call avevar(daty,npty,avey,Sy)  !get average and variance
      v=rnrx+rnry-2.0

      s=Sqrt((rnrx*Sx+rnry*Sy)/v)  !get s (eqn. 5.2)
      
      temp(1)=Sqrt(1.0/rnrx+1.0/rnry)     
      tp=(muxmy)-(avex-avey)/(s*temp(1))  !get tp (eqn 5.2)
      
      temp(2)=rgammln((v-1)/2)-rgammln(v/2)
      temp(3)=1.0/Sqrt(Pi*v)
      
C     calculate probability from eqn 5.3
      buxuy=Exp(temp(2)+log(temp(3)*(1.0+tp*tp/v)**(-(v+1)/2.0)))
     
      return
      end 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotbins(npt,xdata,mean,std,Eerr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,nmax,nplot,i
      parameter(nmax=500 000,nplot=1000)
      real xdata(npt),std,bdatax(nmax),bdatay(nmax),bmax,xb1,xb2,yb1,
     .   yb2,fc,datamin,datamax,mean,plotx(nplot),ploty(nplot),rgauss,
     .   area,Eerr
      
      nbin=real(npt)/20.0
      datamin=mean-4.0*std
      datamax=mean+4.0*std
      call findrange(npt,xdata,datamin,datamax)
      call bindata(nbin,npt,xdata,bdatax,bdatay,datamin,datamax,
     .      bmax)  !bin the data.

      call pgpage() !prepare for multiple pages
      call windowsetup(xb1,xb2,yb1,yb2) !setup a square plot
      call pgvport(xb1,xb2,yb1,yb2) !make a square plotting surface
      fc=2.*(datamax-datamin)/real(nbin)
c      fc=0.
c         write(6,*) "fc:",fc
      call pgwindow(datamin+fc,datamax+fc,0.0,bmax+0.1*bmax)!set window scale
      call pgslw(5)
      call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)!boxes and ticks
      call pglabel("Measurements (mag)","# of Values","") !label axes     
      call pgbin(nbin,bdatax,bdatay,.false.)
      
      area=0.
      do 11 i=1,nbin
         area=area+bdatay(i)*fc/2.0
 11   continue
      
      do 10 i=1,nplot
         plotx(i)=(datamax-datamin)*real(i)/real(nplot)+datamin
         ploty(i)=area*rgauss(std,mean,plotx(i))
 10   continue
 
      call pgsci(2)
      call pgline(nplot,plotx,ploty)
      call pgsci(1)
      
      do 12 i=1,nplot
         plotx(i)=(datamax-datamin)*real(i)/real(nplot)+datamin
         ploty(i)=area*rgauss(Eerr,mean,plotx(i))
 12   continue
 
      call pgsci(3)
      call pgline(nplot,plotx,ploty)
      call pgsci(1)
 
   
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function rgauss(sig,mu,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real sig,mu,x,Pi,temp(2)
      Pi=3.14159265
      
      temp(1)=sig*sqrt(2.0*Pi)
c      temp(1)=1.0
      temp(2)=(x-mu)*(x-mu)/(2*sig*sig)
      rgauss=exp(-temp(2))/temp(1)
      
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

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

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
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
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
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gasdev(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
 1       v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END

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
      real function integrate(npt,x,y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     integrates using Simpson rule.
      implicit none
      integer npt,i,etest
      real x,y(npt)
      
C     start integral
      integrate=y(1)
      do 10 i=2,npt-1
         etest=mod(i,2) !check if value is odd or even
         integrate=integrate+y(i)*(-2.0*real(etest)+4.0)
c         write(*,*) i,etest,y(i),integrate
 10   continue
      integrate=x*(integrate+y(npt))/3.0
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION rgammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     From Numerical Recipes
      real xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      real ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146e0,-86.50532032941677e0,
     * 24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2,
     * -.5395239384953e-5,2.5066282746310005e0/
      x=xx
      y=x
      tmp=x+5.5e0
      tmp=(x+0.5e0)*log(tmp)-tmp
      ser=1.000000000190015e0
      do 11 j=1,6
         y=y+1.e0
         ser=ser+cof(j)/y
 11   continue
      rgammln=tmp+log(stp*ser/x)
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(npt,pdata,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,i
      parameter(nmax=1000)
      real pdata(nmax),datamin,datamax
      
C     initialize datamin and datamax
      datamin=pdata(1)
      datamax=pdata(1)
      
      do 10 i=2,npt
         datamin=min(pdata(i),datamin)
         datamax=max(pdata(i),datamax)
 10   continue
      
      return
      end