C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fabphot(naxes,tix,tiy,fitsdata,itime,gain,rnoise,
     .   fabmag,faberr,sky,a,flat,dpix,nstack,nbegin,nfl)
      implicit none
      integer naxes(2),tix(2),tiy(2),xmax,ymax,npt,i,j,k,nfitm,nfit,
     .   niter,ma,ifun,nn,nstack,nbegin,nsimp,nflsimp,nsky
      parameter(xmax=600,ymax=600,nfitm=200,nsimp=500)
      integer ia(nfitm)
      real fitsdata(xmax,ymax),x1,x2,y1,y2,chi2old,diff,ratio,
     .   pts(xmax*ymax),opx(xmax*ymax),opy(xmax*ymax),avg1,a(nfitm),
     .   x(xmax*ymax),y(xmax*ymax),z(xmax*ymax),alamda,chisq,avg2,nfl,
     .   covar(nfitm,nfitm),alpha(nfitm,nfitm),sig(xmax*ymax),tmp1,tmp2,
     .   itime,gain,fabmag,sky,flux,xhigh,xlow,rnoise,faberr,
     .   opx2(xmax*ymax),opy2(xmax*ymax),flat(xmax,ymax),
     .   dpix(xmax,ymax),simplesum,simpsky,skypix(xmax*ymax),stdsky,
     .   findsky
      save opx,opy,opx2,opy2
      
C     this currently only does flatfield stuff.
c      call fabskyfit(naxes,tix,tiy,fitsdata)


C     go through the data and extract the pixel information for the
C     fabry lense.     
      npt=0
      do 8 j=1,naxes(1)
         do 9 k=1,naxes(2)
            fitsdata(j,k)=fitsdata(j,k)*flat(j,k)+dpix(j,k)
            if((j.ge.tix(1)).and.(j.le.tix(2)).and.
     .           (k.ge.tiy(1)).and.(k.le.tiy(2))) then
               if(flat(j,k).ge.0.6) then
                  npt=npt+1
                  pts(npt)=fitsdata(j,k)
               endif
            endif
 9       continue
 8    continue

C     sort the pixels by intensity      
      call sort(npt,pts)
      do 10 i=1,npt
         x(i)=real(i)
         y(i)=0.
         z(i)=pts(i)
         sig(i)=1.0!sqrt(pts(i))
 10   continue


C     This quick loop summs up the brightest pixels for photometry.
C     e.g. simple photometry..
      simplesum=0.0
      do 11 i=nsimp,npt
         simplesum=simplesum+z(i)
         j=j+1
 11   continue
      nflsimp=j
c      write(0,*) simplesum,z(1),z(npt)
      simpsky=0.0
      do 12 i=1,200
         skypix(i)=z(i)
c         simpsky=simpsky+z(i)
 12   continue
c      simpsky=simpsky/200.0
      nsky=200
      simpsky=findsky(nsky,skypix,stdsky)

C     chop off the 50 brightest pixels (not happy with this!)      
      call rejhilow(npt,x,y,z,0,50)

C     Plotting only first frame 
      if (nbegin.eq.1) then
C     prepare plotting surface      
        call pgvport(0.70,0.95,0.70,0.95)
        x1=0.0
        x2=real(npt)
        y1=0.0
        y2=56384.0*real(nstack)
        call pgsci(0)
C       plot the old intensity distribution (erases it from the screen)
        call plotpoints(npt,opx,opy,1,0,x1,x2,y1,y2)
        call pgline(npt,opx2,opy2)
        call pgsci(1)
        call plotpoints(npt,x,z,1,0,x1,x2,y1,y2)
      endif
      do 20 i=1,npt
         opx(i)=x(i)
         opy(i)=z(i)
 20   continue
 
C     compute averages of the high and low pixels for parameter
C     estimation
      avg1=0.
      do 30 i=1,300
         avg1=avg1+pts(i)
 30   continue
      avg1=avg1/real(300)
      
      avg2=0.
      do 31 i=npt-150,npt
         avg2=avg2+pts(i)
 31   continue 
      avg2=avg2/real(150)

C     Here are the estimated fit parameters (not happy with this)      
      a(3)=600.0!real(npt/2)
      a(1)=(avg1+avg2)/2.0
      a(2)=abs(avg2-avg1)/2.0
      a(4)=50.0

C     set up fitting variables
      nfit=4
      do 35 i=1,nfit
         ia(i)=1
 35   continue
 
      alamda=-1.0
      niter=0
      ma=nfit
      ifun=7
     
C     x - input fitting data
C     y - not used
C     z - input data to fit to
C     sig - error in z
C     npt number of data points
C     a - the fitted parameters
C     ia - which parameters to fit
C     ma - number of fitted parameters
C     covar - work space
C     alpha - work space
C     nfitm - maximum number of allowed fitted values
C     chisq - Chi Squared
C     alamda - check for convergance
C     ifun - fitting function

C     first call to minization routine
      call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,chisq,
     .        alamda,ifun)
      
C     now we iterate minization until convergence
 100  niter=niter+1

c      write(6,*) "niter:",niter,alamda,chisq,(a(i),i=1,nfit)
      chi2old=chisq
      call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,chisq,
     .     alamda,ifun)
      diff =abs(chi2old-chisq)
      ratio=diff/chisq
C     if(stfrm.eq.2) write(6,*) niter,alamda,diff,ratio,chisq
      
      if ((chi2old .gt. chisq .and. ratio .lt. 1.0E-04 .and.
     .     diff .lt. 1.0E-03) .or. (niter .ge. 400) .or.
     .     (alamda .lt. 1.0E-08) .or. (alamda .gt.1.0E08))
     .     then
         alamda=0
         call mrqmin(x,y,z,sig,npt,a,ia,ma,covar,alpha,nfitm,
     .        chisq,alamda,ifun)
      else
         goto 100
      end if

 101  continue


C     plotting for first frame only
      if (nbegin.eq.1)then
C     plot the intensity distribution and the fit      
        do 40 i=1,npt
            y(i)=a(1)+a(2)*tanh((x(i)-a(3))/a(4))
 40     continue
        call pgsci(2)
        call pgline(npt,x,y)
        call pgsci(1)
        do 45 i=1,npt
            opx2(i)=x(i)
            opy2(i)=y(i)
 45     continue
      endif
      
C     okay.. take fitted solution and integrate to estimate star flux
C     and sky flux
C     now we can integate
      xlow=npt-150!a(3)
      xhigh=npt
      tmp1=log(cosh((-a(3)+xlow)/a(4)))
      tmp2=log(cosh((a(3)-xhigh)/a(4)))
      flux=-a(1)*xlow+a(1)*xhigh-a(2)*a(4)*tmp1+a(2)*a(4)*tmp2
      nfl=xhigh-xlow
      
      xlow=0.
      xhigh=200.
      tmp1=log(cosh((-a(3)+xlow)/a(4)))
      tmp2=log(cosh((a(3)-xhigh)/a(4)))
      sky=-a(1)*xlow+a(1)*xhigh-a(2)*a(4)*tmp1+a(2)*a(4)*tmp2
      sky=sky/(xhigh-xlow)
      
c      write(0,*) simplesum,flux
c      write(0,*) simpsky,sky
c      flux=flux-nfl*sky
      flux=simplesum-nflsimp*simpsky

C     The 4 is because the fabry is binned 2x2
      fabmag=25.0-2.5*log10(4.0*flux/itime*gain)
      faberr=1.09/(flux*gain)*sqrt(flux*gain+
     .   4.0*nfl*((rnoise*gain)**2+sky*gain))
c      write(6,*) fabmag,faberr,sky,(a(i),i=1,nfit)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fabskyfit(naxes,tix,tiy,fitsdata)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer naxes(2),tix(2),tiy(2),xmax,ymax,i,j,k,nfitmax,nfit,ma,
     .   mp,np,npt,niter
      parameter(xmax=600,ymax=600,nfitmax=5)
      real fitsdata(xmax,ymax),cvalue,a(nfitmax),u(xmax*ymax,nfitmax),
     .   v(nfitmax,nfitmax),w(nfitmax),zsky,sky,pts(xmax*ymax),
     .   x(xmax*ymax),y(xmax*ymax),z(xmax*ymax),sig(xmax*ymax),chisq

      nfit=nfitmax

      npt=0
      do 8 j=1,naxes(1)
         do 9 k=1,naxes(2)
            if((j.ge.tix(1)).and.(j.le.tix(2)).and.
     .           (k.ge.tiy(1)).and.(k.le.tiy(2))) then
c               if((fitsdata(j,k).ge.datamin)
c     .              .and.(fitsdata(j,k).le.datamax)) then
                  npt=npt+1
                  pts(npt)=fitsdata(j,k)
c               endif
            endif
 9       continue
 8    continue
      
C     sort all the pixel values
      call sort(npt,pts)
C     get the 300th lowest values
      cvalue=pts(300)
      
      npt=0
      do 80 j=1,naxes(1)
         do 90 k=1,naxes(2)
            if((j.ge.tix(1)).and.(j.le.tix(2)).and.
     .           (k.ge.tiy(1)).and.(k.le.tiy(2))) then
               if(fitsdata(j,k).le.cvalue)then
                  npt=npt+1
                  x(npt)=real(j)
                  y(npt)=real(k)
                  z(npt)=fitsdata(j,k)
                  sig(npt)=1.0
               endif
            endif
 90      continue
 80   continue

      ma=nfit
      mp=xmax*ymax
      np=nfit
c      write(6,*) "npt,ma",npt,ma

      niter=5
      zsky=0
      do 52 j=1,niter
         
         call svdfit(x,y,z,sig,npt,a,ma,u,v,w,mp,np,chisq)      
         
         do 50 i=1,npt
            sky=a(1)
            do 51 k=2,ma,2
               sky=sky+a(k)*real(x(i))**(k/2)+
     .              a(k+1)*real(y(i))**(k/2)
 51         continue
            if(j.eq.niter) zsky=zsky+sky
            sig(i)=abs(z(i)-sky)
            if(sig(i).lt.1.0)sig(i)=1.0
 50      continue
         
 52   continue
      zsky=zsky/npt


c      write(6,*) (a(i),i=1,ma)

      do 20 i=1,naxes(1)
         do 21 j=1,naxes(2)
            if((i.ge.tix(1)).and.(i.le.tix(2)).and.
     .           (j.ge.tiy(1)).and.(j.le.tiy(2))) then
               sky=a(1)
               do 22 k=2,ma,2 
                  sky=sky+a(k)*real(i)**(k/2)+
     .                 a(k+1)*real(j)**(k/2)
 22            continue
               write(14,501) i,j,zsky,fitsdata(i,j)
 501  format(2(I4,1X),2(F8.2,1X))
               fitsdata(i,j)=fitsdata(i,j)-sky+zsky
            endif
 21      continue
 20   continue
 
      return
      end
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rejhilow(npt,x,y,z,nlow,nhigh)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Jason Rowe (2005)
      implicit none
      integer npt,nlow,nhigh,ntemp,nmax,i
      parameter(nmax=120000)
      real x(npt),y(npt),z(npt),xt(nmax),yt(nmax),zt(nmax)

      call sort(npt,y,x,z)

      ntemp=0

      do 10 i=1+nlow,npt-nhigh
         ntemp=ntemp+1
         xt(ntemp)=x(i)
         yt(ntemp)=y(i)
         zt(ntemp)=z(i)
 10   continue
      
      do 11 i=1,ntemp
         x(i)=xt(i)
         y(i)=yt(i)
         z(i)=zt(i)
 11   continue
      npt=ntemp

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotpoints(npt,x,y,sym,fy,x1,x2,y1,y2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will plot x,y data values with the symbol defined
C     by sym.  If fy=1 the y-axis is flipped, if fy=0, then it is left
C     alone.
C     if y1=y2 then the y-axis scale is automaticaly calculated
C     same applies for x1,x2
      implicit none
      integer npt,sym,i,fy,ntest,nr
      real x(npt),y(npt),minx,maxx,miny,maxy,x1,x2,y1,y2,px(4),py(4)

C     set these to extremely large and extremely small values.
      minx= 99.9e30
      miny= 99.9e30
      maxx=-99.9e30
      maxy=-99.9e30

C     find the boundaries of the data set
      do 10 i=1,npt
         nr=0
         nr=nr+ntest(x(i))
         nr=nr+ntest(y(i))
         if(nr.eq.0) then
            minx=min(x(i),minx)
            miny=min(y(i),miny)
            maxx=max(x(i),maxx)
            maxy=max(y(i),maxy)
         endif
 10   continue
      
C     Define the boundaries with the data values and add 10%
c      write(6,*) "x1x2:",x1,x2
c      write(6,*) "y1y2:",y1,y2
      if(x1.eq.x2) then
         x1=minx-0.03*(maxx-minx)
         x2=maxx+0.03*(maxx-minx)
      endif
      if(y1.eq.y2) then
         y1=miny-0.03*(maxy-miny)
         y2=maxy+0.03*(maxy-miny)
      endif
c      write(6,*) "x1x2:",x1,x2
c      write(6,*) "y1y2:",y1,y2

C     do we flip the y-axis?
      if(fy.eq.1) then
         call pgwindow(x1,x2,y2,y1)
      else
         call pgwindow(x1,x2,y1,y2)
      endif

      px(1)=x1
      py(1)=y1
      px(2)=x1
      py(2)=y2
      px(3)=x2
      py(3)=y2
      px(4)=x2
      py(4)=y1
c      call pgsci(0)
c      call pgpoly(4,px,py)
c      call pgsci(1)

C     set up the axis and tick marks.
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
C     plots the points!
      call pgpt(npt,x,y,sym)

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
      
