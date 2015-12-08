C Contents
C
C removespline
C removesky
C fitskyspline
C fitskydep

C23456789012345678901234567890123456789012345678901234567890123456789012
C        1         2         3         4         5         6         7

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine removespline(npt,time,mag,sky,itime,mfield,gain,sp,
     .   etimeW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,nph,bins,bsize,i,ntest,nt,sp,steps,panx,pany,
     .   plot,nwbin,jflag
      parameter(nmax=600000)
      real time(npt),mag(npt),sky(npt),gain,x1,x2,xt(nmax),yt(nmax),
     .   yerr(nmax),bx(nmax),by(nmax),itime(npt),etimeW,tmag(nmax),
     .   axx,bxx,cxx,tol,xmin,berr(nmax),fa,fb,fc,oby,nby,f,fmin,brent,
     .   mfield(npt)
      character cans
      common /funct2/ nwbin,bins,bx,by
      
C     sp is for special treatment of hd209458      
                  
      nph=0
      do 60 i=1,npt
         if((itime(i).eq.etimeW).or.(etimeW.eq.0.)) then
c         if(sky(i).gt.10000.0)goto 60
c         if(time(i).gt.2175.25) goto 60
         if(mfield(i).lt.20000.0) goto 60
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C     For HD209458 only
CCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(sp.eq.1) then
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
            elseif(sp.eq.2) then
               if(mfield(i).lt.20000.0) goto 60
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
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END of HD209458 data
CCCCCCCCCCCCCCCCCCCCCCCCCCC
            nph=nph+1
            xt(nph)=sky(i)
            yt(nph)=10**((mag(i)-25.0)/(-2.5))/gain
            yerr(nph)=sqrt(10**((mag(i)-25.0)/(-2.5))/gain)
c            yerr(nph)=sqrt(yt(nph))
c            write(6,*) nph,xt(nph),yt(nph),yerr(nph)
         endif       
 60   continue
c      read(5,*)

c      call sigclip(nph,xt,yt,yerr,3.0)

      write(6,*) "binning data",nph
      if(nph.lt.10) return
C     bins data into equal numbers (so size is dynamic)
      bsize=500
      if(nph/bsize.lt.5) bsize=nph/5
      call bind(nph,xt,yt,yerr,bins,bsize,bx,by,berr)
c      call bindt2(nph,xt,yt,yerr,10.0,3,3)
c      do 61 i=1,nph
c         bx(i)=xt(i)
c         by(i)=yt(i)
c         berr(i)=yerr(i)
c 61   continue
c      bins=nph
      nph=0
      do 1 i=1,bins
         nt=0
         nt=nt+ntest(bx(i))
         nt=nt+ntest(by(i))
         if(nt.eq.0) then
            nph=nph+1
            bx(nph)=bx(i)
            by(nph)=by(i)
            berr(nph)=10.0
         endif
 1    continue
      bins=nph
      

      write(6,*) "Number of bins:",bins

      call skydetrend(npt,mag,sky,itime,gain,bins,bx,by,etimeW)
c      call skyfitnremove(npt,mag,sky,itime,gain,bins,bx,by,berr,etimeW)
            
      return   
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine skyfitnremove(npt,mag,sky,itime,gain,bins,bx,by,berr,
     .   etimeW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,nmax,ma,npc,nfit,i,j,nplot
      parameter(nmax=600000,nfit=5,nplot=1000)
      integer ia(nfit)
      real mag(npt),sky(npt),itime(npt),gain,bx(bins),by(bins),etimeW,
     .   berr(bins),a(nfit),chisq,covar(nfit,nfit),yy,ty,x,
     .   xmin,xmax,ymin,ymax,px(nplot),py(nplot),bxoff,byoff,b(5)


      do 4 i=1,bins
         bx(i)=log(bx(i))
 4    continue
      
      call avgrm(bins,bx,bxoff)
      call avgrm(bins,by,byoff)      
      
      do 5 i=1,nfit
         ia(i)=1
 5    continue
      
      xmin=bx(1)
      xmax=bx(1)
      ymin=by(1)
      ymax=by(1)
      do 10 i=2,bins
         xmin=min(xmin,bx(i))
         xmax=max(xmax,bx(i))
         ymin=min(ymin,by(i))
         ymax=max(ymax,by(i))
 10   continue
      
c      ma=nfit
c      npc=nfit
c      call lfit(bx,by,berr,bins,a,ia,ma,covar,npc,chisq)
c      write(6,*) (a(i),i=1,nfit),chisq
      
c      do 20 i=1,nplot
c         px(i)=xmin+(xmax-xmin)*real(i)/real(nplot)
c         py(i)=0.
c         do 21 j=1,nfit
c            py(i)=py(i)+a(j)*px(i)**(j-1)
c 21      continue
c 20   continue
      
c      call pgwindow(xmin,xmax,ymin,ymax)
c      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
c      call pglabel("Sky","Flux","")
c      call pgpt(bins,bx,by,17)
c      call pgsci(2)
c      call pgline(nplot,px,py)
c      call pgsci(1)
      
      b(1)=2.1502
      b(2)=-0.0603608
      b(3)=-1.97747
      b(4)=6.68149
      b(5)=1.28078e-20
      
      do 30 i=1,npt
         x=sky(i)/4500.0
         yy=b(1)*exp(x*b(2))+b(3)*10**((0.868*b(4)+0.142)*
     .      ((x*b(5))**(1/b(4))))
         yy=-(yy-1.02)*21000.0
         ty=25.0-2.5*log10(yy*gain)
         mag(i)=mag(i)-ty
c         write(6,*) mag(i),ty
 30   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine skydetrend(npt,mag,sky,itime,gain,bins,bx,by,etimeW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,bmax,i
      parameter(bmax=3000)
      real bx(bins),by(bins),y2(bmax),mean,xx,yy,yp1,yp2,
     .     mag(npt),ty,sky(npt),gain,itime(npt),etimeW

      if(bins.gt.bmax) then
         write(6,*) "Increase bmax in skydetrend to:",bins
         pause
      endif
      yp1=1.0e30
      yp2=1.0e30
      call spline(bx,by,bins,yp1,yp2,y2)

      xx=bx(1)
      call splint(bx,by,y2,bins,xx,mean)
      write(6,*) "mean:",xx,mean

      do 10 i=1,npt
         if((itime(i).eq.etimeW).or.(etimeW.eq.0.))then
            xx=sky(i)
            call splint(bx,by,y2,bins,xx,yy)
            ty=25.0-2.5*log10(yy*gain)
c            write(6,*) mag(i),ty,itime(i)
            mag(i)=mag(i)-ty
c            write(6,*) "yy:",yy
CF            ty=10**((mag(i)-25.0)/(-2.5))/gain
CF            ty=ty-yy+mean
c            write(6,*) "ty,yy:",ty,yy
CF            mag(i)=25.0-2.5*log10(ty*gain)
CF            if(ty*gain.lt.0) mag(i)=99.9
         endif
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine removesky(npt,time,mag,sky,itime,nfitc,nfit,at,
     .   aa,gain,etimeW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfitc,nfit,nmax,nfitmax,i,j,nfit2,k
      parameter(nmax=600000,nfitmax=20)
      integer ia(nfitmax)
      real time(npt),mag(npt),at(nfitc),aa(nmax,nfitmax),xx(nmax),
     .   yy(nmax),ans(nfitmax),xerr(nmax),af(nfitmax),gain,itime(nmax),
     .   covar(nfitmax,nfitmax),chisq,x1,x2,y1,y2,ans2(nfitmax,nfitmax),
     .   fli,flo,sky(npt),merr(nmax),avg,adev,sdev,sigma2,skew,curt,
     .   etimeW,ztime

C     order of polynomial to fit.
      nfit2=4
        
      do 1 i=1,nfitc
         xerr(i)=1.0
 1    continue
       
      do 2 i=1,nfit2
         ia(i)=1
 2    continue
      
C     we fit the data for each parameter in aa.      

C     this is done only for convience in the naming convention.
C     pointers would probably be better!
      ztime=at(1)
      do 5 i=1,nfitc    
         xx(i)=at(i)-ztime
 5    continue

c      call pgpage()
C     nfit is the number of parameters we find the sky-mag relation
      do 10 i=1,nfit
C     nfitc is the number of times the relationship was fit.
C     ie. the number of stray light occurances in the data.
         do 15 j=1,nfitc
            yy(j)=aa(j,i)
 15      continue 
         call rejhilow(nfitc,xx,yy,xerr,2,2)
C        make sure we have enough points to fit.
         if(nfitc.lt.nfit2+1) then
            write(6,*) "not enough points",nfitc
            return
         endif
         call lfit(xx,yy,xerr,nfitc,ans,ia,nfit2,covar,nfitmax,chisq)
         x1=0.
         x2=0.
         y1=0
         y2=0
         call pgpage()
         call plotpoints(nfitc,xx,yy,17,0,x1,x2,y1,y2)
         call plotfit(nfit2,ans,x1,x2)
         write(6,*) (ans(j),j=1,nfit2)
         do 16 j=1,nfit2
            ans2(i,j)=ans(j)
 16      continue              
 10   continue
  
      call moment(mag,npt,avg,adev,sdev,sigma2,skew,curt)
      avg=10**((avg-25.0)/(-2.5))/gain

C     for sky removal:
C     sky is in x co-ordinate
C     flux is in y co-ordinate (for fitting purposes) 
C     for parameter estimation
C     x is time co-ordinate, y is the fitted parameter.
C     
C     loop over all data points
      do 20 i=1,npt
C     a different relationship for each exposure time.
         if(itime(i).eq.etimeW) then
C     find interpolated value of fitted parameters
C     I do not care about first fitted parameter of the polynomial.
CC
            do 25 j=2,nfit
               af(j)=0.
               do 26 k=1,nfit2
                  af(j)=af(j)+ans2(j,k)*(time(i)-ztime)**(k-1)
 26            continue
 25         continue
c         write(6,*) (af(j),j=2,nfit),time(i)
            fli=10**((mag(i)-25.0)/(-2.5))/gain
            flo=0.
CC
            do 27 j=2,nfit
               flo=flo+af(j)*sky(i)**(j-1)
 27         continue
C        remove flux here (do I divide or subtract?)
            fli=fli-flo!+avg
            mag(i)=25.0-2.5*log10(fli*gain)
         endif
 20   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitskyspline(npt,sky,flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is not complete!  May never be...
      implicit none
      integer npt,nmax,bins,bsize,i
      parameter(nmax=600000)
      real sky(npt),flux(npt),bx(nmax),by(nmax),berr(nmax),ys(nmax),
     .     x1,x2,y1,y2

C     not using weights.  At least not yet!
      do 10 i=1,npt
         berr(i)=1.0
 10   continue

C     first lets bin the data and remove errendt data points
C     bsize determines howmany datapoints go into each bin
      if(npt.gt.5*15) then
         bsize=npt/15
         call bind(npt,sky,flux,berr,bins,bsize,bx,by,berr)
C         write(6,*) "Bins:",bins
      endif

      call getspline(bins,bx,by,ys)
      x1=0.
      x2=0.
      y1=-10.0
      y2= 10.0
      call plotpoints(bins,bx,ys,17,0,x1,x2,y1,y2)
      write(6,*) (ys(i),i=1,bins)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitskydep(npt,sky,flux,nfit,ans)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      This routine fits a polynomial to the sky/flux relation and 
C      returns the result
      implicit none
      integer npt,nfit,nmax,bins,bsize,nfitmax,ma,i
      parameter(nmax=600000,nfitmax=20)
      integer ia(nfitmax)
      real sky(npt),flux(npt),ans(nfit),bx(nmax),by(nmax),berr(nmax),
     .     covar(nfitmax,nfitmax),chisq,x1,x2,y1,y2
      
      do 5 i=1,nfit
         ia(i)=1
 5    continue
      
C     first lets bin the data and remove errendt data points
C     bsize determines howmany datapoints go into each bin
      bsize=npt/15
      call bind(npt,sky,flux,berr,bins,bsize,bx,by,berr)
C      write(6,*) "Bins:",bins

      do 10 i=1,bins
c         write(6,*) "berr",i,berr(i)
         berr(i)=1.0
 10   continue

C     now we find the polynomial
      ma=nfit
      call lfit(bx,by,berr,bins,ans,ia,ma,covar,nfitmax,chisq)

      x1=0.
      x2=0.
      y1=0.
      y2=0.

C     this was only to test and make sure everything was working
c      call pgpage()
c      call plotpoints(bins,bx,by,17,0,x1,x2,y1,y2)
c      call plotfit(nfit,ans,x1,x2)
c      write(6,*) (ans(i),i=1,nfit)

      return
      end

