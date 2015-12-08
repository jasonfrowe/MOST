CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program mostfit4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     New and improved MOST photometric fixer-upper
      implicit none
      integer nmax,nunit,npt,iargc,nfit,nfitmax,nfitsky,nfitx,nfity,i,j,
     .  nrej,nrejold,niter,nitermax
      parameter(nmax=600000,nfitmax=15)
      integer nd(nmax)
      real px(nmax),py(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),sky(nmax),
     .  xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),tboard(nmax),
     .  mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),itime(nmax),
     .  skystd(nmax),pars(nfitmax),bchisq,chisq,a(nfitmax,nfitmax),
     .  b(nfitmax),diff(nmax),weight(nmax),sigclip
      character*80 filename,outname,cline
c      data pars /12.486,-3.22678d-6,7.00431d-10,-1.71116d-14/
      
      if(iargc().lt.2) goto 901 !get command line arguments. 
      call getarg(1,filename)
      call getarg(2,outname)
      
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,
     .   fy,fxy,tboard,mfield,ftotal,aflux,nap,itime,skystd)
      close(nunit)


C     Default Fitting parameters      
      nfitsky=4 !Sky-magnitude order
      nfitx=3   !x-pixel-magnitude order
      nfity=3   !y-pixel-magnitude order
      sigclip=3.0d0 !outlier detection
      nitermax=10 !number of iterations for robust fitting

C     Read in optional commandline parameters     
      if(iargc().ge.3)then !Sky-magnitude order
        call getarg(3,cline)
        read(cline,*) nfitsky
      endif

      if(iargc().ge.4)then !x-pixel-magnitude order
        call getarg(4,cline)
        read(cline,*) nfitx
      endif
      
      if(iargc().ge.5)then !y-pixel-magnitude order
        call getarg(5,cline)
        read(cline,*) nfity
      endif
      
      if(iargc().ge.6)then !outlier detection
        call getarg(6,cline)
        read(cline,*) sigclip
      endif
      
      if(iargc().ge.7)then !iterations for Robust fitting
        call getarg(7,cline)
        read(cline,*) nitermax
      endif
      
      
      nfit=nfitsky+nfitx+nfity
      if(nfit.gt.nfitmax) then
        write(0,*) "Seg Fault Catch!"
        write(0,*) "nfitmax must be at least :",nfit
      endif
      
      call lfit(nfit,nfitsky,nfitx,nfity,pars,npt,time,mag,merr,sky,xc,
     .  yc,nfitmax,nmax,a,b)  !linear fitting routine 
     
      nrejold=0 !count number of rejected datum
      nrej=1    !initialization of variable
      niter=0 !iteration counter
      do 11 while(nrejold.ne.nrej)
        niter=niter+1 !count number of interations
        nrejold=nrej !save previous count of rejected points
        do 10 i=1,npt !make a copy of photometric errors
            weight(i)=merr(i)
 10     continue     
        call rejectoutlier(nfit,nfitsky,nfitx,pars,npt,time,mag,weight,
     .      sky,xc,yc,diff,sigclip,nrej)
        call lfit(nfit,nfitsky,nfitx,nfity,pars,npt,time,mag,weight,sky,
     .      xc,yc,nfitmax,nmax,a,b)
        if(niter.ge.nitermax) nrejold=nrej !stop if too many iters
 11   enddo
      
c      bchisq=chisq(nfit,nfitsky,pars,npt,time,mag,merr,sky)
      
      call pgopen('?') !open PGPlot device
      call pgask(.false.) !don't ask for new page.. just do it.
C     Sets the size of the plot.  if you have a big monitor increase 6
      call PGPAP ( 6.0 ,1.0) 
      call pgsubp(2,2) !four panels on screen
      call plotsol(nmax,px,py,npt,nd,time,mag,sky,xc,yc,nfit,nfitsky,
     .  nfitx,pars) !plot the solution
      call pgclos()
      
      open(unit=nunit,file=outname) !export detrended data
      call exportdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,fy,fxy,
     .  tboard,mfield,ftotal,aflux,nap,itime,skystd,nfit,nfitsky,nfitx,
     .  pars)
      close(nunit)
      
      goto 999
 901  write(0,500) "Usage: mostfit4 <phot.dat> <outphot.dat> [nfitsky] [
     .nfitx] [nfity] [sigcut] [nitermax]"
 500  format(A87)
      write(0,*) "Where <phot.dat> is a MOST photometry data file"
      write(0,*) "      <outphot.dat> is output name"
      write(0,*) "      [nfitsky] Polynomial Order for Sky-relation fit"
      write(0,*) "      [nfitx] Polynomial Order for X-pixel-relation"
      write(0,*) "      [nfity] Polynomial Order for Y-pixel-relation"
      write(0,*) "      [nitermax] Number of iterations for Robust fitti
     .ng"
      write(0,*) " "
      write(0,*) "Parameters in <> are required, [] are optional" 
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rejectoutlier(nfit,nfitsky,nfitx,pars,npt,time,mag,
     .  merr,sky,xc,yc,diff,sigclip,nrej)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nfitsky,nfitx,npt,i,nrej
      double precision pars(nfit),time(npt),mag(npt),merr(npt),sky(npt),
     .  xc(npt),yc(npt),diff(npt),ans,ave,var,std,func,sigclip
      
      do 10 i=1,npt 
        ans=func(nfit,nfitsky,nfitx,pars,sky(i),xc(i),yc(i)) !model
        diff(i)=mag(i)-ans  !data-model
 10   continue     
      
      call avevar(diff,npt,ave,var)  !get standard deviation 
      std=sqrt(var)
c      write(6,*) "diff: ",ave,var,std
c      read(5,*)
 
      nrej=0
      do 11 i=1,npt
        if(abs(diff(i)-ave).gt.sigclip*std) then !rejection of outliers
            merr(i)=0.0d0
            nrej=nrej+1
        endif
c        write(6,*) i,diff(i),merr(i)
c        read(5,*)
 11   continue
      write(0,*) "nrej:",nrej
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lfit(nfit,nfitsky,nfitx,nfity,pars,npt,time,mag,merr,
     .  sky,xc,yc,nfitmax,nmax,a,b)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,nfitmax,nmax,i,j,k,nfitsky,nfitx,nfity
      double precision pars(nfit),time(npt),mag(npt),merr(npt),
     .  sky(npt),a(nfitmax,nfitmax),b(nfitmax),afunc,xc(npt),yc(npt)
     
C     Solve linear best fit matrix a*x=b     
     
C     set up "a" matrix
      do 10 k=1,nfit
        do 11 j=1,nfit
            a(k,j)=0.0d0
            do 12 i=1,npt
                if(merr(i).gt.0.0d0) a(k,j)=a(k,j)+afunc(nfitsky,nfitx,
     .              nfity,k,sky(i),xc(i),yc(i),merr(i))*afunc(nfitsky,
     .              nfitx,nfity,j,sky(i),xc(i),yc(i),merr(i))
 12         continue
 11     continue
 10   continue
 
c      do 13 i=1,nfit
c        write(0,*) (a(i,j),j=1,nfit)
c 13   continue

C     set up 'b' matrix
      do 15 k=1,nfit
        b(k)=0.0d0
        do 14 i=1,npt
            if(merr(i).gt.0.0d0) b(k)=b(k)+mag(i)*afunc(nfitsky,nfitx,
     .          nfity,k,sky(i),xc(i),yc(i),merr(i))/merr(i)
 14     continue
 15   continue
c      write(0,*) (b(i),i=1,nfit)

C     solve for 'x'
      call gaussj(a,nfit,nfitmax,b,1,1)
c      write(0,*) (b(i),i=1,nfit)
c      read(5,*)
      
      do 16 i=1,nfit
        pars(i)=b(i)
 16   continue
     
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function afunc(nfitsky,nfitx,nfity,n,sky,xc,yc,
     .  merr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,nfitsky,i,nfitx,nfity,n2,n3
      double precision sky,merr,xc,yc

C     This defines the fitting function      
      n2=n-nfitsky
      n3=n-nfitsky-nfitx
      if(n.le.nfitsky) then
        afunc=sky**(n-1)/merr
      elseif(n2.eq.1)then
        if(nfitx.gt.0) afunc=xc/merr
      elseif(n2.eq.2)then
        if(nfitx.gt.1) afunc=int(xc+0.5d0)/merr
      elseif((n2.ge.3).and.(n2.le.nfitx))then
        do 10 i=3,nfitx
            if(nfitx.gt.i-1) afunc=((xc-int(xc+0.5))**(i-1))/merr
 10     continue
      elseif(n3.eq.1)then
        if(nfity.gt.0) afunc=yc/merr
      elseif(n3.eq.2)then
        if(nfity.gt.1) afunc=int(yc+0.5d0)/merr
      elseif((n3.ge.3).and.(n3.le.nfity))then
        do 11 i=3,nfity
            if(nfity.gt.i-1) afunc=((yc-int(yc+0.5))**(i-1))/merr
 11     continue
      endif
c      if(n.gt.nfitsky+nfitx)then
c        write(6,*) n,afunc,yc/merr
c        read(5,*)
c      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotsol(nmax,px,py,npt,nd,time,mag,sky,xc,yc,nfit,
     .  nfitsky,nfitx,pars)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,npt,nfit,i,nplot,nfitsky,nfitx,nd(npt)
      parameter(nplot=1000)
      real px(nmax),py(nmax),bb(4),dsky,dxc,dyc
      double precision time(npt),mag(npt),sky(npt),pars(npt),func,
     .  tmp(2),xc(npt),funcx,medsky,medxc,yc(npt),medyc
     

C     Find medians 
      call rqsort(npt,sky,nd)
      medsky=sky(nd(npt/2))
      call rqsort(npt,xc,nd)
      medxc=xc(nd(npt/2))
      call rqsort(npt,yc,nd)
      medyc=yc(nd(npt/2))
     
    
      bb(1)=real(sky(1))
      bb(2)=real(sky(1))
      bb(3)=real(mag(1))
      bb(4)=real(mag(1))
      do 10 i=1,npt
        px(i)=real(sky(i))
        py(i)=real(mag(i))
        bb(1)=min(px(i),bb(1))
        bb(2)=max(px(i),bb(2))
        bb(3)=min(py(i),bb(3))
        bb(4)=max(py(i),bb(4))
 10   continue

      call pgpage()  
      call pgslw(1)
      call pgsch(2.0)
      call pgvport(0.2,1.0,0.2,0.9)
      call pgwindow(bb(1),bb(2),bb(4),bb(3))
      call pglabel("log(Sky (ADU))","Instrumental Magnitude","")
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pgpt(npt,px,py,-1)
      
      dsky=(bb(2)-bb(1))/real(nplot)
      px(1)=bb(1)
      tmp(1)=dble(px(1))
      tmp(2)=func(nfit,nfitsky,nfitx,pars,tmp(1),medxc,medyc)
      py(1)=real(tmp(2))
      do 11 i=2,nplot
        px(i)=px(i-1)+dsky
        tmp(1)=dble(px(i))
        tmp(2)=func(nfit,nfitsky,nfitx,pars,tmp(1),medxc,medyc)
        py(i)=real(tmp(2))
c        write(0,*) px(i),py(i),medxc
c        read(5,*)
 11   continue
 
      call pgsci(2)
      call pgline(nplot,px,py)
      call pgsci(1)
      
CCCCCCCCCCCCCCCCCCCCCCCC
C     Plot x-co-ordinate
CCCCCCCCCCCCCCCCCCCCCCCC
      bb(1)=real(xc(1))
      bb(2)=real(xc(1))
      bb(3)=real(mag(1))
      bb(4)=real(mag(1))
      do 12 i=1,npt
        px(i)=real(xc(i))
        py(i)=real(mag(i))
        bb(1)=min(px(i),bb(1))
        bb(2)=max(px(i),bb(2))
        bb(3)=min(py(i),bb(3))
        bb(4)=max(py(i),bb(4))
 12   continue

      call pgpage()  
      call pgslw(1)
      call pgsch(2.0)
      call pgvport(0.2,1.0,0.2,0.9)
      call pgwindow(bb(1),bb(2),bb(4),bb(3))
      call pglabel("X-position (Pixel)","Instrumental Magnitude","")
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pgpt(npt,px,py,-1)
      
      dxc=(bb(2)-bb(1))/real(nplot)
      px(1)=bb(1)
      tmp(1)=dble(px(1))
      tmp(2)=func(nfit,nfitsky,nfitx,pars,medsky,tmp(1),medyc)
      py(1)=real(tmp(2))
      do 13 i=2,nplot
        px(i)=px(i-1)+dxc
        tmp(1)=dble(px(i))
        tmp(2)=func(nfit,nfitsky,nfitx,pars,medsky,tmp(1),medyc)
        py(i)=real(tmp(2))
c        write(0,*) px(i),py(i)
c        read(5,*)
 13   continue
 
      call pgsci(2)
      call pgline(nplot,px,py)
      call pgsci(1)      

CCCCCCCCCCCCCCCCCCCCCCCC
C     Plot y-co-ordinate
CCCCCCCCCCCCCCCCCCCCCCCC
      bb(1)=real(yc(1))
      bb(2)=real(yc(1))
      bb(3)=real(mag(1))
      bb(4)=real(mag(1))
      do 14 i=1,npt
        px(i)=real(yc(i))
        py(i)=real(mag(i))
        bb(1)=min(px(i),bb(1))
        bb(2)=max(px(i),bb(2))
        bb(3)=min(py(i),bb(3))
        bb(4)=max(py(i),bb(4))
 14   continue

      call pgpage()  
      call pgslw(1)
      call pgsch(2.0)
      call pgvport(0.2,1.0,0.2,0.9)
      call pgwindow(bb(1),bb(2),bb(4),bb(3))
      call pglabel("Y-position (Pixel)","Instrumental Magnitude","")
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pgpt(npt,px,py,-1)
      
      dyc=(bb(2)-bb(1))/real(nplot)
      px(1)=bb(1)
      tmp(1)=dble(px(1))
      tmp(2)=func(nfit,nfitsky,nfitx,pars,medsky,medxc,tmp(1))
      py(1)=real(tmp(2))
      do 15 i=2,nplot
        px(i)=px(i-1)+dyc
        tmp(1)=dble(px(i))
        tmp(2)=func(nfit,nfitsky,nfitx,pars,medsky,medxc,tmp(1))
        py(i)=real(tmp(2))
c        write(0,*) px(i),py(i)
c        read(5,*)
 15   continue
 
      call pgsci(2)
      call pgline(nplot,px,py)
      call pgsci(1)    
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function func(nfit,nfitsky,nfitx,pars,sky,xc,yc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,nfitsky,nfitx
      double precision pars(nfit),sky,xc,yc
      
      func=0.0d0
      do 10 i=1,nfitsky
        func=func+pars(i)*sky**(i-1)
c        write(6,*) "func:",pars(i),func!pars(i)*sky**(i-1)
10    continue

      do 11 i=1,nfit-nfitsky
        if(i.eq.1) then
            func=func+pars(nfitsky+i)*xc
c            write(6,*) "func:",pars(i),func!pars(i)*sky**(i-1)
        elseif(i.eq.2) then 
            func=func+pars(nfitsky+i)*int(xc+0.5d0)
        else
            func=func+pars(nfitsky+i)*(xc-int(xc+0.5d0))**(i-1)
        endif
11    continue

      do 12 i=1,nfit-nfitsky-nfitx
        if(i.eq.1) then
            func=func+pars(nfitsky+nfitx+i)*yc
c            write(6,*) "func:",pars(i),func!pars(i)*sky**(i-1)
        elseif(i.eq.2) then 
            func=func+pars(nfitsky+nfitx+i)*int(yc+0.5d0)
        else
            func=func+pars(nfitsky+nfitx+i)*(yc-int(yc+0.5d0))**(i-1)
        endif
12    continue

      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function chisq(nfit,nfitsky,pars,npt,time,mag,
     .  merr,sky)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,i,nfitsky
      double precision time(npt),mag(npt),merr(npt),sky(npt),func,chi1,
     .  pars(nfit)
      
      chisq=0.0d0
      
      do 10 i=1,npt
        chi1=(mag(i)-func(nfit,nfitsky,pars,sky))/merr(i)
        chisq=chisq+chi1*chi1
 10   continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,fy,fxy,
     .  tboard,mfield,ftotal,aflux,nap,itime,skystd,nfit,nfitsky,nfitx,
     .  pars)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,nfit,i,nfitsky,nfitx
      double precision time(npt),mag(npt),merr(npt),sky(npt),xc(npt),
     .  yc(npt),fx(npt),fy(npt),fxy(npt),tboard(npt),mfield(npt),
     .  ftotal(npt),aflux(npt),nap(npt),itime(npt),skystd(npt),
     .  pars(nfit),func,ans
     
      do 10 i=1,npt
        ans=func(nfit,nfitsky,nfitx,pars,sky(i),xc(i),yc(i))
        mag(i)=mag(i)-ans
         write(nunit,510) time(i),mag(i),merr(i),
     .      10.0**(sky(i))*itime(i),
     .      xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     .      ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
 10   continue
 
 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
     
      return
      end
      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata(nunit,npt,time,mag,merr,sky,xc,yc,fx,
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
      integer npt,i,ntest,rn
      double precision time(nmax),mag(nmax),merr(nmax),
     .     sky(nmax),xc(nmax),yc(nmax),fx(nmax),fy(nmax),fxy(nmax),
     .     tboard(nmax),mfield(nmax),ftotal(nmax),aflux(nmax),nap(nmax),
     .     itime(nmax),skystd(nmax)

      i=1

   9  continue
      if(i.gt.1) write(6,*) "error reading",i
  10  read(nunit,*,err=9,end=20) time(i),mag(i),merr(i),sky(i),
     . xc(i),yc(i),fx(i),fy(i),fxy(i),tboard(i),mfield(i),
     . ftotal(i),aflux(i),nap(i),itime(i),skystd(i)
        sky(i)=log10(sky(i)/itime(i))
        merr(i)=abs(merr(i))
        if(merr(i).le.0.0) goto 10
        rn=ntest(mag(i))
        if(rn.eq.1) goto 10
        rn=ntest(sky(i))
        if(rn.eq.1) goto 10
C        
C       This is where the outlier rejection routine will be!! 
C

c         if((yc(i).gt.11.0).or.(yc(i).lt.9.0)) goto 10

c        if((fx(i).lt.1.2d0).and.(fy(i).lt.1.3d0).and.
c     .      (fxy(i).gt.-0.03d0).and.(fxy(i).lt.0.06d0).and.
c     .      (tboard(i).gt.292.0d0).and.(mfield(i).gt.20000.0d0).and.
c     .      (yc(i).gt.39.0d0).and.(yc(i).lt.45.0d0).and.
c     .      (xc(i).gt.10.0d0).and.(xc(i).lt.15.0d0))then
            i=i+1
c        endif
      goto 10
 20   continue

 510  format(F13.8,1X,2(F9.6,1X),F9.2,1X,2(F6.2,1X),2(F5.3,1X),F6.3,
     .     1X,F7.3,1X,F9.2,2(1X,F10.2),1X,F8.3,1X,F6.2,1X,F8.4)
      npt=i-1
      write(6,*) "Points read: ",npt

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function ntest(rr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this routine is used to test for NANs and other nasty numbers.
      implicit none
      double precision rr,thi,tlow

C     this defines the numerical boundary of anything that might be 
C     exciting.

      thi=  99.9d30
      tlow=-99.9d30

      if((rr.gt.tlow).and.(rr.lt.thi))then
         ntest=0
      else
         ntest=1
      endif

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER m,mp,n,np,NMAX 
      REAL*8 a(np,np),b(np,mp) 
      PARAMETER (NMAX=1000) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), 
     *     ipiv(NMAX) 
      REAL*8 big,dum,pivinv 
      do 11 j=1,n 
         ipiv(j)=0 
 11   enddo 
      do 22 i=1,n 
         big=0. 
         do 13 j=1,n 
            if(ipiv(j).ne.1)then 
               do 12 k=1,n 
                  if (ipiv(k).eq.0) then 
                     if (abs(a(j,k)).ge.big)then 
                        big=abs(a(j,k)) 
                        irow=j 
                        icol=k 
                     endif 
                  else if (ipiv(k).gt.1) then 
                     write(0,*) "sing:",k
                     stop!pause 'singular matrix in gaussj' 
                  endif 
 12            enddo 
            endif 
 13      enddo 
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then 
            do 14 l=1,n 
               dum=a(irow,l) 
               a(irow,l)=a(icol,l) 
               a(icol,l)=dum 
 14         enddo
            do 15 l=1,m 
               dum=b(irow,l) 
               b(irow,l)=b(icol,l) 
               b(icol,l)=dum 
 15         enddo 
         endif 
         indxr(i)=irow 
         indxc(i)=icol 
         if (a(icol,icol).eq.0.) then
            write(6,*) "sing2:",icol
            stop!pause 'singular matrix in gaussj'
         endif 
         pivinv=1./a(icol,icol) 
         a(icol,icol)=1. 
         do 16 l=1,n 
            a(icol,l)=a(icol,l)*pivinv 
 16      enddo 
         do 17 l=1,m 
            b(icol,l)=b(icol,l)*pivinv 
 17      enddo 
         do 21 ll=1,n 
            if(ll.ne.icol)then 
               dum=a(ll,icol) 
               a(ll,icol)=0. 
               do 18 l=1,n 
                  a(ll,l)=a(ll,l)-a(icol,l)*dum 
 18            enddo
               do 19 l=1,m 
                  b(ll,l)=b(ll,l)-b(icol,l)*dum 
 19            enddo 
            endif 
 21      enddo
 22   enddo 
      do 24 l=n,1,-1 
         if(indxr(l).ne.indxc(l))then 
            do 23 k=1,n 
               dum=a(k,indxr(l)) 
               a(k,indxr(l))=a(k,indxc(l)) 
               a(k,indxc(l))=dum 
 23         enddo 
         endif 
 24   enddo 
      return 
      END 

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
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