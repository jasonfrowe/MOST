C     Contents
C
C     plotdata
C     plotph
C     opengraphics
C     plotpoints
C
C     789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotdata(npt,time,mag,merr,x1,x2,panx,pany,sym)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will plot the data (ie. time vs mag)
C     integer npt: Total number of data points
C     real time(npt): The x-axis data points
C     real mag(npt): The y-axis data points
C     real merr(npt): Error associated with mag()
C     real x1,x2: range of data to plot
C     panx,pany: PGPLOT panel to plot data in

      implicit none
      integer npt,panx,pany,i,sym
      real time(npt),mag(npt),merr(npt),x1,x2,mmin,mmax,tmin,tmax

      call pgpanl(panx,pany)

      mmin= 99.9e30
      mmax=-99.9e30
      tmin= 99.9e30
      tmax=-99.9e30

      do 10 i=1,npt
         if(mag(i).gt.mmax) mmax=mag(i)
         if(mag(i).lt.mmin) mmin=mag(i)
         if(time(i).gt.tmax) tmax=time(i)
         if(time(i).lt.tmin) tmin=time(i)
 10   continue

      if(x1.ne.-1.0) then
         tmin=x1
      endif
      if(x2.ne.-1.0) then
         tmax=x2
      endif

c      mmax=0.1
c      mmin=-0.1
      call pgwindow(tmin,tmax,mmax,mmin)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel("JD","magnitude"," ")
      call pgpt(npt,time,mag,sym)
c      call pgerrb(6,npt,time,mag,merr,1.0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(npt,time,mag,period,panx,pany,sym)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     period - fixed period of the data
C     panx,pany - pgplot panel to plot light curve
      implicit none
      integer nmax
      parameter(nmax=600000)

      integer npt,panx,pany,i,sym
      real time(npt),mag(npt),period,phase(nmax),mmax,mmin,ax1,ax2
      real temp,stdev,std,sigcut,mean


      sigcut=10.0
      mmin= 99.9e30
      mmax=-99.9e30

      mmax=-10.0e10
      mmin= 10.0e10

      mean=0.
      do 5 i=1,npt
         mean=mean+mag(i)
 5    continue
      mean=mean/real(npt)

      std=stdev(npt,mag,mean)
      call pgpanl(panx,pany)
      call phasept(npt,time,mag,phase,period)
      do 10 i=1,npt
         if((mag(i).gt.mmax).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
         if((mag(i).lt.mmin).and.
     .        (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
 10   continue

      ax2=mmin-0.10*(mmax-mmin)
      ax1=mmax+0.10*(mmax-mmin)
      call pgwindow(0.0,1.0,ax1,ax2)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pglabel("Phase","magnitude"," ")
      call pgpt(npt,phase,mag,sym)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotfit(nfit,ans,x1,x2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     plots a polynomial to the already open terminal
C     format is  y = ans(1) + ans(2)*x + ans(3)*x*x...
      implicit none
      integer nfit,nplot,i,j
      parameter(nplot=1000)
      real ans(nfit),x1,x2,xp(nplot),yp(nplot)

      do 10 i=1,nplot
         xp(i)=x1+(x2-x1)*real(i)/real(nplot)
         yp(i)=0.
         do 11 j=1,nfit
            yp(i)=yp(i)+ans(j)*xp(i)**(j-1)
 11      continue
 10   continue

      call pgline(nplot,xp,yp)
         

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine opengraphics(terminal)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine will open the graphics terminal for plotting
C     The character input string is not working very well.
      implicit none
      character*(*) terminal

      call pgopen(terminal)
      call pgpage()
      
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
      real x(npt),y(npt),minx,maxx,miny,maxy,x1,x2,y1,y2

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

C     set up the axis and tick marks.
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
C     plots the points!
      call pgpt(npt,x,y,sym)

      return
      end
